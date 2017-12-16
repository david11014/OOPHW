/*******************************************************
NCKU Department of Mechanical engineering OOP Homework 3
Write by david1104
github: https://github.com/david11014
********************************************************/
#include "HW3_2.h"
#include <memory>

template<class KeyType>
inline Queue<KeyType>::Queue(int MaxQueueSize) : MaxSize(MaxQueueSize + 1), front(0), rear(0)
{
	queue = new KeyType[MaxSize];
};

template<class KeyType>
Queue<KeyType>::Queue(Queue<KeyType>& q) : MaxSize(q.MaxSize), front(q.front), rear(q.rear)
{
	queue = new KeyType[MaxSize];

	for (int i = 0; i < MaxSize; i++)
		queue[i] = q.queue[i];

}

template<class KeyType>
Queue<KeyType>::~Queue()
{
	delete[] queue;
}

template<class KeyType>
Queue<KeyType> Queue<KeyType>::operator=(Queue<KeyType>& q)
{
	return Queue<KeyType>(q);
}

template<class KeyType>
bool Queue<KeyType>::IsFull()
{
	if ((rear + 1) % (MaxSize) == front)
		return true;
	else
		return false;
}

template<class KeyType>
bool Queue<KeyType>::IsEmpty()
{
	if (front == rear)
		return true;
	else
		return false;
}

template<class KeyType>
void Queue<KeyType>::Add(const KeyType & item)
{
	if (!IsFull())
	{
		queue[rear] = item;
		rear = (rear + 1) % (MaxSize);
	}
	else
		QueueFull();

}

template<class KeyType>
KeyType * Queue<KeyType>::Delete()
{
	if (!IsEmpty())
	{
		KeyType* t = &(queue[front]);

		front = (front + 1) % (MaxSize);

		return t;
	}
	else
	{
		QueueEmpty();
		return 0;
	}
		


}
