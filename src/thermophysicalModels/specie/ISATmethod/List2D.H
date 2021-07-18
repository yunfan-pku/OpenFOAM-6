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
    int timeSteps_;
    int size_;
    int NtimeSteps_;
    std::vector<Node> list_;
    //std::vector<Node*> listtail;
    int start_, end_;

public:
    inline List2D(int n) : timeSteps_(0), size_(0), NtimeSteps_(n), list_(n), start_(0), end_(n - 1)
    {
        for (int i = 0;i < NtimeSteps_;i++)
        {
            list_[i].next = &list_[i];
            list_[i].prev = &list_[i];
            list_[i].data = nullptr;
            list_[i].timetag = -1;
        }
    }
    inline const int& NtimeSteps()const
    {
        return NtimeSteps_;
    }

    inline std::vector<Node>& list()
    {
        return list_;
    }

    ~List2D()
    {
        print_simple();

        for (int i = 0;i < NtimeSteps_;i++)
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
        temp->timetag = timeSteps_;
        temp->prev = list_[start_].prev;
        temp->next = &list_[start_];
        list_[start_].prev->next = temp;
        list_[start_].prev = temp;
        size_++;
        return temp;
    }
    inline T* pop()
    {
        if (size_ == 0)
            return nullptr;
        int i = end_;
        for (;list_[i].prev == &list_[i];i = (i + NtimeSteps_ - 1) % NtimeSteps_);
        Node* pnode = list_[i].prev;
        T* temp = pnode->data;
        pnode->prev->next = pnode->next;
        pnode->next->prev = pnode->prev;
        size_--;
        delete pnode;
        return temp;
    }
    inline void renew(Node* x)
    {
        if (x->timetag == timeSteps_)
            return;
        x->prev->next = x->next;
        x->next->prev = x->prev;
        x->prev = list_[start_].prev;
        x->next = &list_[start_];
        list_[start_].prev->next = x;
        list_[start_].prev = x;
        x->timetag = timeSteps_;
    }
    inline void newTimeStep()
    {
        timeSteps_++;
        int newend_ = (end_ + NtimeSteps_ - 1) % NtimeSteps_;
        if (list_[end_].next != list_[end_].prev)
        {
            list_[newend_].prev->next = list_[end_].next;
            list_[end_].next->prev = list_[newend_].prev;
            list_[newend_].prev = list_[end_].prev;
            list_[end_].prev->next = &list_[newend_];
            list_[end_].prev = &list_[end_];
            list_[end_].next = &list_[end_];
        }
        start_ = (start_ + NtimeSteps_ - 1) % NtimeSteps_;
        end_ = newend_;
    }
    void print()
    {
        int j;
        std::cout << "size_:" << size_ << std::endl;
        for (int i = 0;i < NtimeSteps_;i++)
        {
            std::cout << i << " :";
            j=(i + start_) % NtimeSteps_;
            for (Node* p = &list_[j];p->next != &list_[j];p = p->next)
            {
                std::cout << "(" << p->next->data->value()[0] << "," << p->next->data->value()[1] << ")" << "->";
            }
            std::cout << std::endl;
        }

    }
    inline void print_simple()
    {
        int sum;
        int j;
        std::cout << "size_:" << size_ << std::endl;
        for (int i = 0;i < NtimeSteps_;i++)
        {
            sum = 0;
            std::cout << i << " :";
            j=(i + start_) % NtimeSteps_;
            for (Node* p = &list_[j];p->next != &list_[j];p = p->next)
            {
                sum++;
                //std::cout << "("<<p->next->data->value()[0]<<","<< p->next->data->value()[1]<<")"<< "->";
            }
            std::cout << sum << std::endl;
        }

    }

    inline void leaflist(T* retlist[])
    {
        int iter = 0;
        int j;
        for (int i = 0;i < NtimeSteps_;i++)
        {
            j = (i + start_) % NtimeSteps_;

            for (Node* p = &list_[j];p->next != &list_[j];p = p->next)
            {
                retlist[iter++] = p->next->data;
            }
        }
        //std::cout << iter  << std::endl;

    }
};

#endif