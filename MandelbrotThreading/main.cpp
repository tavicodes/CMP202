// Mandelbrot set example
// Adam Sampson <a.sampson@abertay.ac.uk>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>

#include "structs.h"
#include "farm.h"
#include "mandelbrotTask.h"

#include <thread>
#include <mutex>

#include <iomanip>
#include <time.h>
#include <string>
#include <numeric>

#include <amp.h>
#include <amp_math.h>

// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::complex;
using std::cout;
using std::cin;
using std::endl;
using std::ofstream;

using namespace concurrency;

// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_clock;

// The mutex to lock the image data 
std::mutex imageMutex;

// The image data.
// Each pixel is represented as 0xRRGGBB.
uint32_t image[HEIGHT][WIDTH];

/* CPU CODE */

// Render the Mandelbrot set into the image array.
// The parameters specify the region on the complex plane to plot.
void compute_mandelbrotCPU(mandelbrotRegion reg, imageSizeCPU imSize)
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
				std::unique_lock<std::mutex> lock(imageMutex);
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = 0x000000; // black
			}
			else
			{
				std::unique_lock<std::mutex> lock(imageMutex);
				// z escaped within less than reg.maxIterations
				// iterations. This point isn't in the set.
				image[y][x] = 0xFFFFFF; // white
			}
		}
	}
}

// Set up CPU threading for the Mandelbrot set to be computed
// The parameters specify the region on the complex plane to plot and whether to split the strips into rows or columns
void start_mandelbrotCPU(mandelbrotRegion reg, int mode)
{
	imageSizeCPU stripSize;
	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	std::vector<std::thread*> threadVector;

	if (mode == 0)
	{
		stripSize.maxHeight = 0;
		int heightDifference = HEIGHT / concurentThreadsSupported;

		for (int i = 0; i < (int)concurentThreadsSupported; i++)
		{
			stripSize.maxHeight += heightDifference;
			threadVector.push_back(new std::thread(compute_mandelbrotCPU, reg, stripSize));
			stripSize.minHeight = stripSize.maxHeight;
		}
	}
	else if (mode == 1)
	{
		stripSize.maxWidth = 0;
		int widthDifference = WIDTH / concurentThreadsSupported;
		for (int i = 0; i < (int)concurentThreadsSupported; i++)
		{
			stripSize.maxWidth += widthDifference;
			threadVector.push_back(new std::thread(compute_mandelbrotCPU, reg, stripSize));
			stripSize.minWidth = stripSize.maxWidth;
		}
	}
	else
	{
		Farm cpuFarm;

		/*uint32_t* pointImRows[HEIGHT];
		for (int i = 0; i < HEIGHT; i++)
		{
			pointImRows[i] = image[i];
		}
		uint32_t** fullPointIm = pointImRows;*/

		int heightDifference = HEIGHT;
		int widthDifference = WIDTH;
		int heightEdit = 0;
		int widthEdit = 0;
		
		bool oddsFound = false;

		while (!oddsFound)
		{
			if (heightDifference % 2 != 1)
			{
				heightDifference = heightDifference / 2;
				heightEdit++;
			}
			else if (widthDifference % 2 != 1)
			{
				widthDifference = widthDifference / 2;
				widthEdit++;
			}
			else
			{
				oddsFound = true;
			}
		}

		cout << "Setting up task farm..." << endl;

		for (int currentHeight = 0; currentHeight < (1 << heightEdit); currentHeight++)
		{
			stripSize.minHeight = currentHeight * heightDifference;
			stripSize.maxHeight = (currentHeight + 1) * heightDifference;

			for (int currentWidth = 0; currentWidth < (1 << widthEdit); currentWidth++)
			{
				stripSize.minWidth = currentWidth * heightDifference;
				stripSize.maxWidth = (currentWidth + 1) * heightDifference;

				cpuFarm.add_task(new mandelbrotTask(reg, stripSize));
				
			}
		}
		
		cout << "Running task farm..." << endl;

		cpuFarm.run();
	}

	if (mode < 2)
	{
		for (int i = 0; i < (int)concurentThreadsSupported; i++)
		{
			threadVector[i]->join();
		}
	}
}

/* GPU CODE */

// using our own Complex number structure and definitions as the Complex type is not available in the Concurrency namespace
struct Complex1 {
	float x;
	float y;
};

Complex1 c_add(Complex1 c1, Complex1 c2) restrict(cpu, amp)	// restrict keyword - able to execute this function on the GPU and CPU
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a + c;
	tmp.y = b + d;
	return tmp;
}// c_add

float c_abs(Complex1 c) restrict(cpu, amp)
{
	return concurrency::fast_math::sqrt(c.x * c.x + c.y * c.y);
}// c_abs

Complex1 c_mul(Complex1 c1, Complex1 c2) restrict(cpu, amp)
{
	Complex1 tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a * c - b * d;
	tmp.y = b * c + a * d;
	return tmp;
}// c_mul

void report_accelerator(const accelerator a)
{
	const std::wstring bs[2] = { L"false", L"true" };
	std::wcout << ": " << a.description << " "
		<< endl << "       device_path                       = " << a.device_path
		<< endl << "       dedicated_memory                  = " << std::setprecision(4) << float(a.dedicated_memory) / (1024.0f * 1024.0f) << " Mb"
		<< endl << "       has_display                       = " << bs[a.has_display]
		<< endl << "       is_debug                          = " << bs[a.is_debug]
		<< endl << "       is_emulated                       = " << bs[a.is_emulated]
		<< endl << "       supports_double_precision         = " << bs[a.supports_double_precision]
		<< endl << "       supports_limited_double_precision = " << bs[a.supports_limited_double_precision]
		<< endl;
}
// List and select the accelerator to use
void list_accelerators()
{
	//get all accelerators available to us and store in a vector so we can extract details
	std::vector<accelerator> accls = accelerator::get_all();

	// iterates over all accelerators and print characteristics
	for (unsigned i = 0; i < accls.size(); i++)
	{
		accelerator a = accls[i];
		report_accelerator(a);
		//if ((a.dedicated_memory > 0) & (a.dedicated_memory < 0.5*(1024.0f * 1024.0f)))
		//accelerator::set_default(a.device_path);
	}

	accelerator::set_default(accls[0].device_path);
	accelerator acc = accelerator(accelerator::default_accelerator);
	std::wcout << " default acc = " << acc.description << endl;

} // list_accelerators

  // query if AMP accelerator exists on hardware
void query_AMP_support()
{
	std::vector<accelerator> accls = accelerator::get_all();
	if (accls.empty())
	{
		cout << "No accelerators found that are compatible with C++ AMP" << std::endl;
	}
	else
	{
		cout << "Accelerators found that are compatible with C++ AMP" << std::endl;
		list_accelerators();
	}
} // query_AMP_support

// set accelerator to use
// 0 is GPU, 3 is CPU
void set_accelerator(int a)
{
	//get all accelerators available to us and store in a vector
	std::vector<accelerator> accls = accelerator::get_all();
	accelerator::set_default(accls[a].device_path);
}

// Render the Mandelbrot set into the image array.
// The parameters specify the region on the complex plane to plot.
void compute_mandelbrotGPU(mandelbrotRegion reg)
{
	const int iterationMax = reg.maxIterations;

	uint32_t* pImage = &(image[0][0]);
	array_view<uint32_t, 2> a(HEIGHT, WIDTH, pImage);
	a.discard_data();

	parallel_for_each(a.extent, [=](concurrency::index<2> idx) restrict(amp) {
		// compute Mandelbrot here i.e. Mandelbrot kernel/shader
		//USE THREAD ID/INDEX TO MAP INTO THE COMPLEX PLANE
		int y = idx[0];
		int x = idx[1];

		// Work out the point in the complex plane that
		// corresponds to this pixel in the output image.
		Complex1 c;
		c.x = reg.left + (x * (reg.right - reg.left) / WIDTH);
		c.y = reg.top + (y * (reg.bottom - reg.top) / HEIGHT);

		// Start off z at (0, 0).
		Complex1 z;
		z.x = 0.0;
		z.y = 0.0;

		// Iterate z = z^2 + c until z moves more than 2 units
		// away from (0, 0), or we've iterated too many times.
		int iterations = 0;
		while (c_abs(z) < 2.0 && iterations < iterationMax)
		{
			z = c_add(c_mul(z, z), c);

			++iterations;
		}

		if (iterations == iterationMax)
		{
			// z didn't escape from the circle.
			// This point is in the Mandelbrot set.
			a[y][x] = 0x000000; // black
		}
		else
		{
			// z escaped within less than reg.maxIterations
			// iterations. This point isn't in the set.
			a[y][x] = 0xFFFFFF; // white
		}
		});
}

/* MAIN BODY CODE */

// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char* filename)
{
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char*)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char*)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}

int main(int argc, char* argv[])
{
	query_AMP_support();

	std::vector<long long> timeCPURow;
	std::vector<long long> timeCPUColumn;
	std::vector<long long> timeAMPGPU;
	std::vector<long long> timeAMPCPU;
	mandelbrotRegion setRegion;

	std::string mandelbrotRegionCheck;

	cout << "Would you like to display a portion of the Mandelbrot set? (Y/N)" << endl;
	cin >> mandelbrotRegionCheck;

	if (mandelbrotRegionCheck == "Y" || mandelbrotRegionCheck == "y")
	{
		cout << "Please enter the co-ordinates of the Mandelbrot region:" << endl << "Left: ";
		cin >> setRegion.left;
		cout << "Right: ";
		cin >> setRegion.right;
		cout << "Top: ";
		cin >> setRegion.top;
		cout << "Bottom: ";
		cin >> setRegion.bottom;
		cout << "Please enter the maximum number of iterations value: ";
		cin >> setRegion.maxIterations;
	}
	else
	{
		setRegion.left = -2.0;
		setRegion.right = 1.0;
		setRegion.top = 1.125;
		setRegion.bottom = -1.125;
	}

	cout << "Running computations..." << endl;

	for (int i = 0; i < 10; i++)
	{
		cout << "CPU (rows) running..." << endl;

		// Start timing, compute full mandelbrot and stop timing
		the_clock::time_point start = the_clock::now();
		start_mandelbrotCPU(setRegion, 0);
		the_clock::time_point end = the_clock::now();

		// Compute the difference between the two times in milliseconds and add to vector
		auto time_taken = duration_cast<microseconds>(end - start).count();
		timeCPURow.push_back(time_taken);

		cout << "Computing the Mandelbrot set with CPU (rows) took " << time_taken << " us." << endl;

		cout << "CPU (columns) running..." << endl;

		// Start timing, compute full mandelbrot and stop timing
		start = the_clock::now();
		//start_mandelbrotCPU(setRegion, 1);
		end = the_clock::now();

		// Compute the difference between the two times in milliseconds and add to vector
		time_taken = duration_cast<microseconds>(end - start).count();
		timeCPUColumn.push_back(time_taken);

		cout << "Computing the Mandelbrot set with CPU (columns) took " << time_taken << " us." << endl;

		write_tga("outputFarm.tga");

		cout << "CPU (task farm) running..." << endl;

		// Start timing, compute full mandelbrot and stop timing
		start = the_clock::now();
		start_mandelbrotCPU(setRegion, 2);
		end = the_clock::now();

		// Compute the difference between the two times in milliseconds and add to vector
		time_taken = duration_cast<microseconds>(end - start).count();
		timeCPUColumn.push_back(time_taken);

		cout << "Computing the Mandelbrot set with task farm (chunks) took " << time_taken << " us." << endl;

		cout << "AMP GPU running..." << endl;

		// Set accelerator to GPU
		set_accelerator(0);

		// Start timing, compute full mandelbrot and stop timing
		start = the_clock::now();
		compute_mandelbrotGPU(setRegion);
		end = the_clock::now();

		// Compute the difference between the two times in milliseconds and add to vector
		time_taken = duration_cast<microseconds>(end - start).count();
		timeAMPGPU.push_back(time_taken);

		cout << "Computing the Mandelbrot set with AMP GPU took " << time_taken << " us." << endl;

		cout << "AMP CPU running..." << endl;

		// Set accelerator to CPU
		set_accelerator(1);

		// Start timing, compute full mandelbrot and stop timing
		start = the_clock::now();
		compute_mandelbrotGPU(setRegion);
		end = the_clock::now();

		// Compute the difference between the two times in milliseconds and add to vector
		time_taken = duration_cast<microseconds>(end - start).count();
		timeAMPCPU.push_back(time_taken);

		cout << "Computing the Mandelbrot set with AMP CPU took " << time_taken << " us." << endl;
	}

	cout << "Writing to timings.txt \n";

	ofstream timingsFile;
	timingsFile.open("timings.txt");
	timingsFile << "CPU Row timings: \n";
	for (int i = 0; i < timeCPURow.size(); i++)
	{
		timingsFile << timeCPURow[i] << std::endl;
	}
	timingsFile << "CPU Column timings: \n";
	for (int i = 0; i < timeCPUColumn.size(); i++)
	{
		timingsFile << timeCPUColumn[i] << std::endl;
	}
	timingsFile << "AMP GPU timings: \n";
	for (int i = 0; i < timeAMPGPU.size(); i++)
	{
		timingsFile << timeAMPGPU[i] << std::endl;
	}
	timingsFile << "AMP CPU timings: \n";
	for (int i = 0; i < timeAMPCPU.size(); i++)
	{
		timingsFile << timeAMPCPU[i] << std::endl;
	}
	timingsFile.close();

	cout << "Writing successful. \n";

	// This zooms in on an interesting bit of detail.
	//compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);
	//write_tga("outputTest.tga");

	return 0;
}
