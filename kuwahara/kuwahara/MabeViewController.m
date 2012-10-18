//
//  MabeViewController.m
//  kuwahara
//
//  Created by Sweeple on 10/18/12.
//  Copyright (c) 2012 mabesoft. All rights reserved.
//

#import "MabeViewController.h"
#import "UIImage+kuwahara.h"

@interface MabeViewController ()

@end

@implementation MabeViewController

- (void)viewDidLoad
{
    [super viewDidLoad];
	// Do any additional setup after loading the view, typically from a nib.
    _filteredImageView_r2.image = [_originalImageView.image imageByApplyingKuwaharaFilterWithRadius:2];
    _filteredImageView_r3.image = [_originalImageView.image imageByApplyingKuwaharaFilterWithRadius:3];
    _filteredImageView_r4.image = [_originalImageView.image imageByApplyingKuwaharaFilterWithRadius:4];
    _filteredImageView_r6.image = [_originalImageView.image imageByApplyingKuwaharaFilterWithRadius:6];
    _filteredImageView_r8.image = [_originalImageView.image imageByApplyingKuwaharaFilterWithRadius:8];
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (void)dealloc {
    [_originalImageView release];
    [_filteredImageView_r2 release];
    [_filteredImageView_r3 release];
    [_filteredImageView_r4 release];
    [_filteredImageView_r6 release];
    [_filteredImageView_r8 release];
    [super dealloc];
}

@end
