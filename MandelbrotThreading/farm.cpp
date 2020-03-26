#include "farm.h"

#include <thread>
#include <vector>

void Farm::add_task(Task* task)
{
	taskMutex.lock();
	taskQueue.push(task);
	taskMutex.unlock();
}

void Farm::run()
{
	int threadCount = std::thread::hardware_concurrency();
	std::vector<std::thread*> threadVector;

	for (int i = 0; i < threadCount; i++)
	{
		threadVector.push_back(new std::thread([&] {
			while (true)
			{
				if (taskQueue.empty())
				{
					return;
				}

				taskMutex.lock();
				Task* currentTask = taskQueue.front();
				taskQueue.pop();
				taskMutex.unlock();

				currentTask->run();
				delete currentTask;
			}
			}));
	}

	for (int i = 0; i < threadCount; i++)
	{
		threadVector[i]->join();
	}
}