#ifndef LIST2D_H
#define LIST2D_H
#include<vector>
#include<iostream>

template<class T>
class List2D
{
public:
    struct Node
    {
        Node* next;
        Node* prev;
        T* data;
        int timetag;
    };
private:
    int timeSteps;
    int size;
    int NtimeSteps;
    std::vector<Node> list_;
    //std::vector<Node*> listtail;
    int start, end;

public:
    inline List2D(int n) : timeSteps(0), size(0), NtimeSteps(n), list_(n), start(0), end(n - 1)
    {
        for (int i = 0;i < NtimeSteps;i++)
        {
            list_[i].next = &list_[i];
            list_[i].prev = &list_[i];
            list_[i].data = nullptr;
            list_[i].timetag = -1;
        }
    }
    ~List2D()
    {
        print_simple();

        for (int i = 0;i < NtimeSteps;i++)
        {
            for (Node* p = list_[i].next->next;p->prev != &list_[i];p = p->next)
            {
                //std::cout << "1! ";
                //std::cout << (p->prev) << " " << &list_[i];
                delete p->prev;
                //std::cout << "2!" << std::endl;
                if (p == &list_[i])
                    break;
            }
        }
        //std::cout << "1!!" << std::endl;
    }
    inline Node* insert(T* x) {
        Node* temp = new Node;
        temp->data = x;
        temp->timetag = timeSteps;
        temp->prev = list_[start].prev;
        temp->next = &list_[start];
        list_[start].prev->next = temp;
        list_[start].prev = temp;
        size++;
        return temp;
    }
    inline T* pop()
    {
        if (size == 0)
            return nullptr;
        int i = end;
        for (;list_[i].prev == &list_[i];i = (i + NtimeSteps - 1) % NtimeSteps);
        Node* pnode = list_[i].prev;
        T* temp = pnode->data;
        pnode->prev->next = pnode->next;
        pnode->next->prev = pnode->prev;
        size--;
        delete pnode;
        return temp;
    }
    inline void renew(Node* x)
    {
        if (x->timetag == timeSteps)
            return;
        x->prev->next = x->next;
        x->next->prev = x->prev;
        x->prev = list_[start].prev;
        x->next = &list_[start];
        list_[start].prev->next = x;
        list_[start].prev = x;
        x->timetag = timeSteps;
    }
    inline void newTimeStep()
    {
        timeSteps++;
        int newend = (end + NtimeSteps - 1) % NtimeSteps;
        if (list_[end].next != list_[end].prev)
        {
            list_[newend].prev->next = list_[end].next;
            list_[end].next->prev = list_[newend].prev;
            list_[newend].prev = list_[end].prev;
            list_[end].prev->next = &list_[newend];
            list_[end].prev = &list_[end];
            list_[end].next = &list_[end];
        }
        start = (start + NtimeSteps - 1) % NtimeSteps;
        end = newend;
    }
    void print()
    {
        std::cout << "size:" << size << std::endl;
        for (int i = 0;i < NtimeSteps;i++)
        {
            std::cout << i << " :";
            for (Node* p = &list_[(i + start) % NtimeSteps];p->next != &list_[(i + start) % NtimeSteps];p = p->next)
            {
                std::cout << "(" << p->next->data->value()[0] << "," << p->next->data->value()[1] << ")" << "->";
            }
            std::cout << std::endl;
        }

    }
    void print_simple()
    {
        int sum;
        std::cout << "size:" << size << std::endl;
        for (int i = 0;i < NtimeSteps;i++)
        {
            sum = 0;
            std::cout << i << " :";
            for (Node* p = &list_[(i + start) % NtimeSteps];p->next != &list_[(i + start) % NtimeSteps];p = p->next)
            {
                sum++;
                //std::cout << "("<<p->next->data->value()[0]<<","<< p->next->data->value()[1]<<")"<< "->";
            }
            std::cout << sum << std::endl;
        }

    }
};

#endif