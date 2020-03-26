#pragma once

// The size of the image to generate.
const int WIDTH = 1920;
const int HEIGHT = 1200;

// Struct for Mandelbrot set region
struct mandelbrotRegion
{
	double left;
	double right;
	double top;
	double bottom;

	// The number of times to iterate before we assume that a point isn't in the
	// Mandelbrot set.
	// (You may need to turn this up if you zoom further into the set.)
	int maxIterations = 500;
};

// Struct for CPU thread image size
struct imageSizeCPU
{
	int minWidth = 0;
	int maxWidth = WIDTH;
	int minHeight = 0;
	int maxHeight = HEIGHT;
};