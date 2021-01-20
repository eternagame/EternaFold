#ifndef __INC_RAND_H__
#define __INC_RAND_H__

#include <cstdlib>
class Die
{
public:
  Die(unsigned int seed)
  {
    srand(seed);
  }

  double operator()()
  {
    return rand()/(RAND_MAX+1.0);
  }
};

template < class T, class RealT >
class Roulette
{
public:
    Roulette(Die &die)
        : t_(), max_(NEG_INF), die_(die)
    { }

    void add(T t, RealT v)
    {
        t_.push_back(std::make_pair(t, v));
        max_ = std::max(max_, v);
    }

    T choose() const
    {
        RealT sum = 0.0;
        std::list<std::pair<T,RealT> > u;
        typename std::list< std::pair<T,RealT> >::const_iterator x;
        for (x=t_.begin(); x!=t_.end(); ++x)
        {
            u.push_back(std::make_pair(x->first, Fast_Exp(x->second - max_)));
            sum += u.back().second;
        }

        Assert(sum>0.0, "sum not > 0");
        RealT r = die_()*sum;
        RealT s = 0.0;
        for (x=u.begin(); x!=u.end(); ++x)
        {
            s += x->second;
            if (r<s) return x->first;
        }
        return u.back().first;
    }

    unsigned int size() const { return t_.size(); }

private:
    std::list<std::pair<T,RealT> > t_;
    RealT max_;
    Die& die_;
};
#endif  //  __INC_RAND_H__
