#ifndef HW3_2_H
#define HW3_2_H

#include <iostream>
#define DefaultQueueSize 100

using namespace std;
template <class KeyType> class Queue;

template <class KeyType>
ostream& operator<<(ostream& os, const Queue<KeyType>& p);

template <class KeyType>
class Queue
{
private:
	int front, rear;
	KeyType* queue;
	int MaxSize;

	void QueueFull() { cout << "Queue is Full!!\n"; };
	void QueueEmpty() { cout << "Queue is Empty!!\n"; };

public:
	Queue(int MaxQueueSize = DefaultQueueSize);
	// Create an empty queu whose maximum size is MaxQueueSize
	Queue(Queue<KeyType>& );
	~Queue();

	Queue<KeyType> operator= (Queue<KeyType>&);

	bool IsFull();
	// if number of elements in the queue is equal to the maximum size of
	// the queue, return TRUE; otherwise, return FALSE
	bool IsEmpty();
	// if number of elements in the queue is equal to 0, return TRUE
	void Add(const KeyType& item);
	// if IsFull(), then QueueFull(); else insert item at rear of the queue
	KeyType* Delete();
	// if IsEmpty(), then QueueEmpty() and return 0;
	// else remove the item at the front of the queue and return a pointer to it
	
	friend ostream& operator<< <>(ostream& os, const Queue& p);
};

template <class KeyType>
ostream& operator<<(ostream& os, const Queue<KeyType>& q)
{	
	for (int i = q.front; i != q.rear ; i = (i + 1) % q.MaxSize)
	{
		os << q.queue[i] << " ";
	}

	return os;
}

#endif // !HW3_2_H