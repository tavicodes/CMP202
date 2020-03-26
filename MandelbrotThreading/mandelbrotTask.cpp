#include "mandelbrotTask.h"

void mandelbrotTask::run()
{
	for (int y = imSize.minHeight; y < imSize.maxHeight; ++y)
	{
		for (int x = imSize.minWidth; x < imSize.maxWidth; ++x)
		{
			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			complex<double> c(reg.left + (x * (reg.right - reg.left) / WIDTH),
				reg.top + (y * (reg.bottom - reg.top) / HEIGHT));

			// Start off z at (0, 0).
			complex<double> z(0.0, 0.0);

			// Iterate z = z^2 + c until z moves more than 2 units
			// away from (0, 0), or we've iterated too many times.
			int iterations = 0;
			while (abs(z) < 2.0 && iterations < reg.maxIterations)
			{
				z = (z * z) + c;

				++iterations;
			}

			if (iterations == reg.maxIterations)
			{
				std::unique_lock<std::mutex> lock(imageTaskMutex);
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = 0x000000; // black
			}
			else
			{
				std::unique_lock<std::mutex> lock(imageTaskMutex);
				// z escaped within less than reg.maxIterations
				// iterations. This point isn't in the set.
				image[y][x] = 0xFFFFFF; // white
			}
		}
	}
}
