//
//  UIImage+kuwahara
//
//  Created on 10/12/12.
//
//  Copyright (c) 2012 by Marco Benedetti
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted
//  provided that the following conditions are met:
//
//  Redistributions of source code must retain the above copyright notice, this list of conditions and the
//  following disclaimer.
//  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
//  following disclaimer in the documentation and/or other materials provided with the distribution.
//  Neither the name of the UIImage+kuwahara software nor the names of its contributors may be used to endorse
//  or promote products derived from this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
//  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
//  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.

#import <UIKit/UIKit.h>

@interface UIImage (kuwahara)
// This method applies a Kuwahara filter with the given radius to an image (self) and returns the filtered image.
//
// Take a look here for a description of how the filter works:
//  http://en.wikipedia.org/wiki/User:DinoVgk
//
// The implementation attempts to be fast by using GCD, Accelerate, and an incremental algorithm.
// This filter does not use the GPU, only the CPU.
-(UIImage*) imageByApplyingKuwaharaFilterWithRadius:(NSUInteger) radius;
@end
