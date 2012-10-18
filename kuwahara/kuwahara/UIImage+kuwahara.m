//
//  UIImage+kuwahara
//
//  Created on 10/12/12.
//
//  Copyright (c) 2012 by Marco Benedetti
//  All rights reserved.

#import "UIImage+kuwahara.h"
#import <Accelerate/Accelerate.h>

// If you comment out the following definition, scalar code will be used instead of the
// vectorized code produced by the "Accelerate" framework (which leverages the NEON
// instruction set of ARM processors); it looks like using NEON speeds up the computation
// of the kawahara filter by at least a factor of 2 (twice as fast or more)
#define USE_ACCELERATE

// If you comment out the following definition, the alpha channel is NOT skipped by the algorithm
// but it takes part in all the computations from the input to the output image; this is often
// a waste of work because images undergoing the kawahara filter have no transparent region.
// Furthermore, the alpha channel plays no role in the criterion used to compute the output value.
// If the alpha channel is ignored, the output image will have exactly the same alpha-values
// of the input one; skipping the alpha channel slighlty improves performance
#define SKIP_ALPHA_CHANNEL

// The image will be divided into #DEGREE_OF_PARALLELISM equally-sized vertical strips
// before being filtered, and each strip will be worked on separately, i.e., "in parallel"
// (GCD is used to parallelize the work, so the actual degree of parallelism depends on
// the number of cores on the target machine, on the computational load, etc.)
#define DEGREE_OF_PARALLELISM (4)

#pragma mark -
#pragma mark === Types used to manipulate the image ===
#pragma mark

typedef struct {
    float a, r, g, b;
} ccmp;

@implementation UIImage (kuwahara)

void _CG_drawImageExtendingBorders(CGContextRef ctx, CGRect rect, CGImageRef img);
void _kuwahara_filterImage(size_t output_width, size_t output_height, ccmp* output, ccmp* input, size_t radius);

-(UIImage*) imageByApplyingKuwaharaFilterWithRadius:(NSUInteger) radius {
    
    // Let's first define the size of the images we'll be working with
    // They'll be larger than self to accomodate for a surrounding border
    // the size of the radius meant to simplify the implementation of the filter
    size_t width  = 2*radius + CGImageGetWidth(self.CGImage);
    size_t height = 2*radius + CGImageGetHeight(self.CGImage);
    
    size_t nPixels = width * height; // number of pixels in the image
    size_t nComponents = nPixels * 4; // number of color components in the image (four components per pixel: A,R,G,B)
    size_t nBytes = nComponents * sizeof(float); // number of bytes required to represent a float-based version of the image
    
    // The context that will contain the input image for the filter
    CGColorSpaceRef colorSpace = CGColorSpaceCreateDeviceRGB();
    CGContextRef input_context = CGBitmapContextCreate (NULL, width, height, 8, 4*width, colorSpace,
#ifdef SKIP_ALPHA_CHANNEL
    // if we decided to skip the alpha channel, we don't need such channel in the context; furthermore, we need to
    // work with the non-premultiplied values of the r,g,b channels, and this settings gives us just that
    kCGImageAlphaNoneSkipFirst
#else
    // if we decided to keep the alpha channel, we are ok with the premultiplied values, so we decide to
    // work with a context in the native iOS image format to speed up things a bit
    kCGImageAlphaPremultipliedFirst|kCGBitmapByteOrder32Little
#endif
);
    CGColorSpaceRelease(colorSpace);
    
    // Then we draw the image (ourself) in the middle of the context, precisely within the border
    CGRect rect = CGRectMake(radius, radius, CGImageGetWidth(self.CGImage), CGImageGetHeight(self.CGImage));
    _CG_drawImageExtendingBorders(input_context, rect, self.CGImage);
        
    // The context is translated from an integer-based to a float-based representation
    // This is necessary because on one hand iOS does not natively support float-based contexts
    // and on the other hand the Accelerate primitives we'll be using work with floats
    void *input_bitmapData_int = CGBitmapContextGetData(input_context);
    ccmp *input_bitmapData_float = (ccmp *)malloc(nBytes);
    vDSP_vfltu8(input_bitmapData_int, 1, (float*)input_bitmapData_float, 1, nComponents);
    
    // We're ready for the actual filtering, performed by the _kuwahara_filterImage function on the float-based image representation
    ccmp *output_bitmapData_float = (ccmp *)malloc(nBytes);
    _kuwahara_filterImage(width, height, input_bitmapData_float, output_bitmapData_float, radius);
#ifdef SKIP_ALPHA_CHANNEL
    // If we decided not to process the alpha channel, the current alpha values
    // in the output image are meaningless; so, we copy those of the input image
    cblas_scopy(nPixels, (float*)input_bitmapData_float, 4, (float*)output_bitmapData_float, 4);
#endif
    free(input_bitmapData_float); // the float input to the filter is no longer used
    
    // The float-based representation just obtained as a result of the filtering operation is converted back to an int-based representation
    // The old integer bitmap data for the input context (not deallocated and no longer in use) is reused to host the int-based representation
    CGContextRef output_context = input_context;
    void *output_bitmapData_int = input_bitmapData_int;
    vDSP_vfixu8((float*)output_bitmapData_float, 1, output_bitmapData_int, 1, nComponents);
    free(output_bitmapData_float); // the float output of the filter is no longer required
    
    // Finally, we remove the artificial border introduced to facilitate the filter and turn the resulting CGImage into a UIImage
    CGImageRef output_image_withBorder = CGBitmapContextCreateImage(output_context);
    CGImageRef output_image_noBorder = CGImageCreateWithImageInRect(output_image_withBorder, rect);
    UIImage *result = [UIImage imageWithCGImage:output_image_noBorder scale:[self scale] orientation:[self imageOrientation]];
    CGImageRelease(output_image_noBorder);
    CGImageRelease(output_image_withBorder);
    CGContextRelease(output_context);
    
    return result;
}

-(UIImage*) imageByApplyingKuwaharaFilterWithRadius4 {
    return [self imageByApplyingKuwaharaFilterWithRadius:3];
}

@end

#pragma mark -
#pragma mark === Addition to CoreGraphics: image drawing with border extension
#pragma mark

void _CG_subimagePattern(void *info, CGContextRef ctx) {
    CGImageRef img = (CGImageRef)info;
    CGRect rect = CGRectMake(0,0,CGImageGetWidth(img),CGImageGetHeight(img));
    CGContextDrawImage(ctx, rect, img);
}

typedef enum {
    SqueezeNot            = 0x00,
    SqueezeTowardMinXEdge = 0x01,
    SqueezeTowardMinYEdge = 0x02,
    SqueezeTowardMaxXEdge = 0x04,
    SqueezeTowardMaxYEdge = 0x08,
    SqueezeTowardTLCorner = SqueezeTowardMinXEdge|SqueezeTowardMinYEdge,
    SqueezeTowardTRCorner = SqueezeTowardMaxXEdge|SqueezeTowardMinYEdge,
    SqueezeTowardBLCorner = SqueezeTowardMinXEdge|SqueezeTowardMaxYEdge,
    SqueezeTowardBRCorner = SqueezeTowardMaxXEdge|SqueezeTowardMaxYEdge
} _CG_SqueezeDirection;

CGRect _CG_squeezeRect(CGRect r, _CG_SqueezeDirection e) {
    if (e&SqueezeTowardMinXEdge || e&SqueezeTowardMaxXEdge) {
        if (e&SqueezeTowardMaxXEdge) r.origin.x += (r.size.width-1.0f);
        r.size.width = 1.0f;
    }
    if (e&SqueezeTowardMinYEdge || e&SqueezeTowardMaxYEdge) {
        if (e&SqueezeTowardMaxYEdge) r.origin.y += (r.size.height-1.0f);
        r.size.height = 1.0f;
    }
    return r;
}

void _CG_drawImageExtendingBorders(CGContextRef ctx, CGRect rect, CGImageRef img) {

    size_t width  = CGBitmapContextGetWidth(ctx);
    size_t height = CGBitmapContextGetHeight(ctx);
    
	static const CGPatternCallbacks callbacks = {0, &_CG_subimagePattern, NULL};
    
    CGSize rSize[3][3] = {
        (CGSize){CGRectGetMinX(rect),       CGRectGetMinY(rect)},  // Top left
        (CGSize){rect.size.width,           CGRectGetMinY(rect)},  // Top center
        (CGSize){width-CGRectGetMaxX(rect), CGRectGetMinY(rect)},
        (CGSize){CGRectGetMinX(rect),       rect.size.height},
        (CGSize){rect.size.width,           rect.size.height},
        (CGSize){width-CGRectGetMaxX(rect), rect.size.height},
        (CGSize){CGRectGetMinX(rect),       height-CGRectGetMaxY(rect)},
        (CGSize){rect.size.width,           height-CGRectGetMaxY(rect)},
        (CGSize){width-CGRectGetMaxX(rect), height-CGRectGetMaxY(rect)}
    };

    _CG_SqueezeDirection squeezeDir[3][3] = {
        SqueezeTowardTLCorner, SqueezeTowardMinYEdge, SqueezeTowardTRCorner,
        SqueezeTowardMinXEdge,            SqueezeNot, SqueezeTowardMaxXEdge,
        SqueezeTowardBLCorner, SqueezeTowardMaxYEdge, SqueezeTowardBRCorner
    };
        
    CGColorSpaceRef patternSpace = CGColorSpaceCreatePattern (NULL);
	CGContextSetFillColorSpace (ctx, patternSpace);
	CGColorSpaceRelease (patternSpace);

    CGRect imgRect = CGRectMake(0,0,CGImageGetWidth(img),CGImageGetHeight(img));
    CGAffineTransform tr = CGAffineTransformMakeTranslation(rect.origin.x, rect.origin.y);
    
	float alpha = 1;
    CGRect targetRect = CGRectZero;
    for (int y=0; y<3; y++, targetRect.origin.y += targetRect.size.height, targetRect.origin.x = 0.0f) {
        for (int x=0; x<3; x++, targetRect.origin.x += targetRect.size.width) {
            targetRect.size = rSize[y][x];
            if (x==1 && y==1) {
                CGContextDrawImage(ctx, targetRect, img);
            } else {
                CGRect sourceRect = _CG_squeezeRect(imgRect,squeezeDir[y][x]);
                CGPatternRef pattern = CGPatternCreate(img, sourceRect, tr, sourceRect.size.width, sourceRect.size.height, kCGPatternTilingNoDistortion, true, &callbacks);
                CGContextSetFillPattern (ctx, pattern, &alpha);
                CGPatternRelease (pattern);
                CGContextFillRect(ctx, targetRect);
            }
        }
    }
}

#pragma mark -
#pragma mark === Data structures used to represent an incremental Kawahara window ===
#pragma mark

typedef struct {
    struct { ccmp mean, sigma; } left;
    struct { ccmp mean, sigma; } right;
} floatVec;

#define EMPTY_VEC {0.0f, 0.0f, 0.0f, 0.0f}
floatVec floatVec0 = { {EMPTY_VEC,EMPTY_VEC}, {EMPTY_VEC,EMPTY_VEC} };

// The kuwahara window, defined not in terms of its absolute position in the image
// but in terms of current sub-totals for the relevant values (mean, sigma) in the
// four sub-windows and in its (1+2*radius) rows.
typedef struct {
    NSUInteger radius; // the radius of the window
    
    floatVec  TOP; // current mean and sigma-square of the values in the top left and right sub-windows
    floatVec  BOT; // current mean and sigma-square of the values in the bottom left and right sub-windows
    floatVec *ROW; // (circular) vector of (1+2*radius) means and sigma-squares for the half-rows on the left and right halves of the window
    NSInteger ROW_midIdx; // the index (within ROW[]) of the row that is currently central
} KuwaharaWindow;

#pragma mark -
#pragma mark === Vector manipulation functions (with/without Accelerate) ===
#pragma mark

CF_INLINE void _vec_accumulate(floatVec *acc, floatVec *data) {
    size_t nFloats = sizeof(floatVec)/sizeof(float); // 16
#ifdef USE_ACCELERATE
    cblas_saxpy(nFloats, 1.0, (float*)data, 1, (float*) acc, 1);
#else
    for (int i=0; i<nFloats; i++)
        ((float*)(acc))[i] += ((float*)(data))[i];
#endif
}

CF_INLINE void _vec_subtract(floatVec *acc, floatVec *data) {
    size_t nFloats = sizeof(floatVec)/sizeof(float); // 16
#ifdef USE_ACCELERATE
    cblas_saxpy(nFloats,-1.0, (float*)data, 1, (float*) acc, 1);
#else
    for (int i=0; i<nFloats; i++)
        ((float*)(acc))[i] -= ((float*)(data))[i];
#endif
}

CF_INLINE ccmp _vec_sum(ccmp *row, size_t rowSize) {
    ccmp result;
#ifdef USE_ACCELERATE
#ifndef SKIP_ALPHA_CHANNEL
    vDSP_sve((float*)row  , 4, &(result.a), rowSize);
#endif
    vDSP_sve((float*)row+1, 4, &(result.r), rowSize);
    vDSP_sve((float*)row+2, 4, &(result.g), rowSize);
    vDSP_sve((float*)row+3, 4, &(result.b), rowSize);
#else
    result.a = result.r = result.g = result.b = 0.0f;
    for (int i=0; i<rowSize; i++) {
#ifndef SKIP_ALPHA_CHANNEL
        result.a += row[i].a;
#endif
        result.r += row[i].r;
        result.g += row[i].g;
        result.b += row[i].b;
    }
#endif
    return result;
}

CF_INLINE ccmp _vec_sumOfSquares(ccmp *row, size_t rowSize) {
    ccmp result;
#ifdef USE_ACCELERATE
#ifndef SKIP_ALPHA_CHANNEL
    vDSP_svesq((float*)row  , 4, &(result.a), rowSize);
#endif
    vDSP_svesq((float*)row+1, 4, &(result.r), rowSize);
    vDSP_svesq((float*)row+2, 4, &(result.g), rowSize);
    vDSP_svesq((float*)row+3, 4, &(result.b), rowSize);
#else
    result.a = result.r = result.g = result.b = 0.0f;
    for (int i=0; i<rowSize; i++) {
#ifndef SKIP_ALPHA_CHANNEL
        result.a += row[i].a*row[i].a;
#endif
        result.r += row[i].r*row[i].r;
        result.g += row[i].g*row[i].g;
        result.b += row[i].b*row[i].b;
    }
#endif
    return result;
}

#pragma mark -
#pragma mark === Kawahara window manipulation functions ===
#pragma mark

KuwaharaWindow *_kuwaharaWindow_allocWithRadius(NSUInteger radius) {
    KuwaharaWindow *w = (KuwaharaWindow *)malloc(sizeof(KuwaharaWindow));
    w->radius = radius;
    w->ROW_midIdx = radius;
    NSUInteger nr = 1+2*radius; ///< number of rows in the window
    w->ROW = (floatVec*)malloc(sizeof(floatVec)*nr);
    return w;
}

void _kuwaharaWindow_release(KuwaharaWindow *w) {
    free(w->ROW);
    free(w);
}

void _kuwaharaWindow_init(KuwaharaWindow *w, ccmp *inboundRow, size_t row_delta) {
    NSUInteger r = w->radius;
    NSUInteger nr = 1+2*(w->radius);
    w->TOP = w->BOT = floatVec0;
    for (NSUInteger i=0; i<nr; i++, inboundRow+=row_delta) {
        (w->ROW)[i] = (floatVec) {
            { _vec_sum(inboundRow,   1+r), _vec_sumOfSquares(inboundRow,   1+r) },
            { _vec_sum(inboundRow+r, 1+r), _vec_sumOfSquares(inboundRow+r, 1+r) }
        };
        if (i<=r) {
            // mean and sigma of the top left and right subwindow are updated
            _vec_accumulate(&(w->TOP),  (w->ROW)+i);
        }
        if (i>=r) {
            // mean and sigma of the botton left and right subwindow are updated
            _vec_accumulate(&(w->BOT),  (w->ROW)+i);
        }
    }
    // Given the way we've processed rows, the mid index is in the middle entry of ROW
    w->ROW_midIdx = w->radius;
}

CF_INLINE NSInteger _kuwaharaWindow_idxOfRow(KuwaharaWindow *w, NSInteger relIdx) {
    NSUInteger nr = 2*(w->radius)+1; // number of rows in the window
    return ((w->ROW_midIdx)+relIdx+nr)%nr; // modular arithmetics to properly index rows in our circular array
}

// This function shifts the kuwahara window w one step/pixel/row down and updates incrementally
// the mean and sigma of the 4 sub-windows plus the sub-totals stored in w->ROW.
// The row of the (1+2*r) inbound pixels that enter the window as a result of the shift
// must be pointed to by the input parameter inboundRow.
// The array ROW is treated as circular, so moving it by one step is just a matter of
// incrementing the reference "mid row index" and to update one single entry (the entry that
// contained the previous top row - which exits the window - will contain the new bottom row
// - which enters the window).
// Similarly, the mean-s and sigma-s are updated by subtracting the values that are leaving the
// window and adding those that are entering it.
// This incremental approach allows for a number of operations-per-pixel proportional to the
// window radius and not to the square of the radius.
void _kuwaharaWindow_shiftDown(KuwaharaWindow *w, ccmp *inboundRow) {
    
    NSUInteger r = w->radius;
    NSUInteger idxOfTopRow = _kuwaharaWindow_idxOfRow(w,-r), idxOfMidRow = _kuwaharaWindow_idxOfRow(w, 0);

    // First, the mean/sigma values of the four subwindows are incrementally updated
    // by subtracting the values related to the outbound row that is exiting the window
    _vec_subtract(&(w->TOP), &((w->ROW)[idxOfTopRow]));
    _vec_subtract(&(w->BOT), &((w->ROW)[idxOfMidRow]));
    
    // The values of the new row of pixels the window is just entering is taken into account
    (w->ROW)[idxOfTopRow] = (floatVec) {
        { _vec_sum(inboundRow,   1+r), _vec_sumOfSquares(inboundRow,   1+r) },
        { _vec_sum(inboundRow+r, 1+r), _vec_sumOfSquares(inboundRow+r, 1+r) },
    };

    // Then, the window is shifted by one row (downward)
    w->ROW_midIdx = _kuwaharaWindow_idxOfRow(w,+1);

    // Now, the indexes of the bottom and mid rows in ROW[] are as follows
    NSUInteger idxOfBotRow = _kuwaharaWindow_idxOfRow(w,+r); idxOfMidRow = _kuwaharaWindow_idxOfRow(w, 0);
    
    // Finally, the mean/sigma values of the four subwindows are incrementally updated
    // by adding the values related to the inbound rows that just entered the window
    _vec_accumulate(&(w->TOP), &((w->ROW)[idxOfMidRow]));
    _vec_accumulate(&(w->BOT), &((w->ROW)[idxOfBotRow]));
}

// Assuming w is a kuwahara window centered at point (x,y) on the input image, this
// function returns the value of the point (x,y) after the kuwahara filter is applied.
// Note that w maintains no reference to the absolute position (x,y) on the image, so
// by "centered at point (x,y) on the input image" we just mean that the 4 means and
// the 4 sigmas in w are indeed the means and sigmas of the four subwindows of a
// kuwahara window centered at (x,y).
ccmp _kuwaharaWindow_currentFilteredValue(KuwaharaWindow *w) {
    float n = (1+w->radius)*(1+w->radius); // Number of pixels in each of the four sub-windows
    float s_min = MAXFLOAT;

#ifdef USE_ACCELERATE
    ccmp m[4] = {w->TOP.left.mean,  w->TOP.right.mean,  w->BOT.left.mean,  w->BOT.right.mean};
    ccmp s[4] = {w->TOP.left.sigma, w->TOP.right.sigma, w->BOT.left.sigma, w->BOT.right.sigma};
    
    ccmp mDivByN[4]; vDSP_vsdiv((float*)m, 1, &n, (float*)&mDivByN, 1, 16);
    ccmp sDivByN[4]; vDSP_vsdiv((float*)s, 1, &n, (float*)&sDivByN, 1, 16);
    ccmp mDivByN2[4]; vDSP_vsq((float*)&mDivByN, 1, (float*)&mDivByN2, 1, 16);
    ccmp mDivByN_minus_mDivByN2[4]; vDSP_vsub((float*)&mDivByN2, 1, (float*)&sDivByN, 1, (float*)&mDivByN_minus_mDivByN2, 1, 16);
    ccmp absOf_mDivByN_minus_mDivByN2[4]; vDSP_vabs ((float*)&mDivByN_minus_mDivByN2, 1, (float*)&absOf_mDivByN_minus_mDivByN2, 1, 16);
    
    // Which one of the four sub-window exhibits the lowest variance?
    // We return the average value of the pixels in that sub-window
    int k_best = 0;
    for (int k = 0; k < 4; k++) {
        float s_sum;
        vDSP_sve(1+(float*)&(absOf_mDivByN_minus_mDivByN2[k]), 1, &s_sum, 3);
        if (s_sum<s_min) {
            s_min = s_sum;
            k_best = k;
        }
    }
    return mDivByN[k_best];
#else
    ccmp result;
    ccmp *m = &(w->TOP.left.mean);
    ccmp *s = &(w->TOP.left.sigma);
    
    // Which one of the four sub-window exhibits the lowest variance?
    // We return the average value of the pixels in that sub-window
    for (int k = 0; k < 4; k++, m+=2, s+=2) {
        ccmp m1 = { 0.0f, (m->r)/n, (m->g)/n, (m->b)/n};
        ccmp s1 = { 0.0f, fabs((s->r)/n-m1.r*m1.r), fabs((s->g)/n-m1.g*m1.g), fabs((s->b)/n-m1.b*m1.b) };
#ifndef SKIP_ALPHA_CHANNEL
        m1.a = (m->a)/n;
        s1.a = fabs((s->a)/n-m1.a*m1.a);
#endif
        float s_sum = s1.r + s1.g + s1.b;
        if (s_sum<s_min) {
            s_min = s_sum;
            result = m1;
        }
    }
    return result;
#endif
}

#pragma mark -
#pragma mark === Kawahara filtering algorithm ===
#pragma mark

void _kuwahara_filterImage(size_t width, size_t height, ccmp* input, ccmp* output, size_t radius) {

    // We'll process the image in parallel dividing it into numberOfStrips vertical strips
    size_t numberOfStrips = DEGREE_OF_PARALLELISM;
    size_t stripWidth = (width-2*radius)/numberOfStrips;
    size_t remainder  = (width-2*radius)%numberOfStrips;

    // This is the group of blocks dispatched to work on the different image strips
    dispatch_group_t group = dispatch_group_create();
    
    // We can't be sure that the image width divided by numberOfStrips is an integer number;
    // in general there will be a remainder (from 0 to numberOfStrips); if R is the reminder,
    // we arrange for the first R strips to be one pixel wider than the remaining ones so we
    // perfectly cover the whole image surface
    stripWidth++;
    for (NSUInteger stripIdx=0, x0=0; stripIdx<numberOfStrips; stripIdx++, x0+=stripWidth) {
        if (stripIdx==remainder)
            stripWidth--;

        dispatch_group_async(group, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^ {
            // Each strip worker allocates an indipendent kuwahara window
            KuwaharaWindow *w = _kuwaharaWindow_allocWithRadius(radius);
            ccmp *inboundRow = NULL;
            // All the columns of the current strip are considered
            for (NSUInteger x=radius+x0, s=0; s<stripWidth; x++, s++) {
                // The window is initially loaded with the pixels at the top of the current column
                inboundRow = input+x-radius;
                _kuwaharaWindow_init(w, inboundRow, width);
                inboundRow += (2*radius)*width;
                // we go through all the rows of the current column and push the window down one row at a time
                for (NSInteger y=radius; y<height-radius; y++, inboundRow += width, _kuwaharaWindow_shiftDown(w, inboundRow)) {
                    // Here the kuwahara window is centered at pixel (x,y) od the image; we ask it to compute
                    // the filtered value of the current pixel and to assign the result to the output buffer
                    output[width*y+x] = _kuwaharaWindow_currentFilteredValue(w);
                }
            }
            _kuwaharaWindow_release(w);
        });
    }
    
    // Let's wait for all the image strips to be processed before returning
    dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
    dispatch_release(group);
}
