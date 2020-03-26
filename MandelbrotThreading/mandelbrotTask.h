#pragma once

#include "task.h"
#include "structs.h"

#include <cstdint>
#include <complex>
#include <mutex>

using std::complex;

extern uint32_t image[HEIGHT][WIDTH];

class mandelbrotTask : public Task
{
public:
	mandelbrotTask(mandelbrotRegion newReg, imageSizeCPU newImSize)
		: reg(newReg), imSize(newImSize)
	{
	}

	void run();
private:
	mandelbrotRegion reg;
	imageSizeCPU imSize;
	std::mutex imageTaskMutex;
};

