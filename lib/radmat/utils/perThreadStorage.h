#ifndef PERTHREADSTORAGE_H_H_GUARD
#define PERTHREADSTORAGE_H_H_GUARD

#include <map>
#include <omp.h>
#include <iostream>

namespace radmat
{

// a scoped lock to unset on exit
struct guardLock
{
    guardLock(omp_lock_t &lock_)
        : lock(&lock_) , own(false)
    {
        acquire();
    }

    ~guardLock(void)
    {
        release();
    }

    void acquire(void)
    {
        omp_set_lock(lock);
        own = true;
    }

    void release(void)
    {
        if(own)
            {
                omp_unset_lock(lock);
                own = false;
            }
    }

    omp_lock_t *lock;
    bool own;
};




// not memory safe if spawning off pointers
template<class T>
struct ompPerThreadMap
{

    typedef typename std::map<int,T>::const_iterator const_iterator; 

    // this guy locks 
    T &operator()(void)
    {
        // scoped lock
        guardLock lock(lock_);

        // try to find the T specific to this thread
        typename std::map<int, T>::iterator it = threadMap.find(omp_get_thread_num());

        // if it exist return it
        if(it != threadMap.end())
            return it->second;

        // insert it and then return the T
        // stl definition -- pair<iterator,bool> insert ( const value_type& x );
        return threadMap.insert(typename std::map<int, T>::value_type(omp_get_thread_num(), T())).first->second;
    }

    // read only !!! 
    T &operator()(const int id)
    {
      typename std::map<int, T>::iterator it = threadMap.find(id);

      // if it exist return it
      if(it == threadMap.end())
        inject_thread(omp_get_thread_num()); 

      it = threadMap.find(id); 

      return it->second;
    }

    void pre_inject_thread(const int id)
    {
      // stl definition -- pair<iterator,bool> insert ( const value_type& x );
      threadMap.insert(typename std::map<int, T>::value_type(id, T()));
    }

    void inject_thread(const int id)
    {
      // scoped lock
      guardLock lock(lock_);

      // stl definition -- pair<iterator,bool> insert ( const value_type& x );
      threadMap.insert(typename std::map<int, T>::value_type(id, T()));
    }

    const_iterator begin(void) const {return threadMap.begin();}
    const_iterator end(void) const {return threadMap.end();}


    void emptyThisThread(void)
    {
      omp_set_lock(&lock_);
      typename std::map<int, T>::iterator it = threadMap.find(omp_get_thread_num());

      if(it != threadMap.end())
        threadMap.erase(it);

      omp_unset_lock(&lock_);
    }

  private:
    std::map<int, T> threadMap;
    omp_lock_t lock_;
};


/*
   test case

   struct boobar
   {

   void doSomething(void)
   {
   threadLocalDouble() = omp_get_thread_num();
   }

   void doSomethingElse(void)
   {
   doSomething();

#pragma omp critical
{
// std::cout << threadLocalDouble() << std::endl;
std::cout << omp_get_thread_num() << std::endl;
}
}

private:
static ompPerThreadMap<double> threadLocalDouble;
};


ompPerThreadMap<double> boobar::threadLocalDouble;

*/

} //radmat

#endif
