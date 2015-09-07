/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2002, 2009, 2010, 2015 Tino Kluge (ttk448 at gmail.com)
 * Copyright (C) 2009, 2010 MathFinance AG (info@mathfinance.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see
 *  <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 */


#ifndef MF_FDM_OPTION_HPP
#define MF_FDM_OPTION_HPP

namespace fin
{

// namespace for enum's
namespace type
{

enum payoff {
    vanilla = 0,
    cash    = 1,
    asset   = 2
};

enum put_call {
    put     = 0,
    call    = 1
};

enum asian {
    fixed_strike    = 0,
    floating_strike = 1
};

} // namespace type


// covers put/call, digital options, knock out barrier and no touch options
struct OptionBarrier {
public:
    type::payoff     pay;       // vanilla, cash/asset or nothing
    type::put_call   putcall;
    double  strike;             // strike
    double  T;                  // time to maturity
    double  barrier_lower;      // possible barriers (ignored, if zero)
    double  barrier_upper;
    double  rebate_lower;       // payment if barrier is reached
    double  rebate_upper;

    // constructor (default)
    OptionBarrier()
    {
        pay     = type::vanilla;
        putcall = type::call;
        barrier_lower = 0.0;
        barrier_upper = 0.0;
        rebate_lower  = 0.0;
        rebate_upper  = 0.0;
    };

    // payoff profile of the option at maturity T
    double payoff(double s) const
    {
        bool has_lower_barrier = (barrier_lower>0.0);
        bool has_upper_barrier = (barrier_upper>0.0);
        double  intrinsic=0.0;
        double  pc = (putcall==type::call) ? 1.0 : -1.0;
        switch(pay) {
        case type::vanilla:
            intrinsic=std::max(pc*(s-strike),0.0);
            break;
        case type::cash:      // digital: cash or nothing
            if((putcall==type::call && s>=strike)
                    || (putcall==type::put && s<=strike) ) {
                intrinsic=1.0;
            } else {
                intrinsic=0.0;
            }
            break;
        case type::asset:      // digital: asset or nothing
            if((putcall==type::call && s>=strike)
                    || (putcall==type::put && s<=strike) ) {
                intrinsic=s;
            } else {
                intrinsic=0.0;
            }
            break;
        default:
            // option payoff not implemented
            assert(false);
        }

        if(has_lower_barrier && s<=barrier_lower)
            intrinsic=rebate_lower;
        if(has_upper_barrier && s>=barrier_upper)
            intrinsic=rebate_upper;
        return intrinsic;
    }
};


struct OptionAsian {
public:
    type::put_call  putcall;
    type::asian     fixed_float;
    double  strike;             // strike, or strike factor
    double  t;                  // current time
    double  T;                  // time of maturity
    std::vector<double> dates;  // averaging schedule
    std::vector<double> spot;   // for past dates, we need to know spot fixes


    // constructor (default)
    OptionAsian()
    {
        putcall = type::call;
        fixed_float = type::fixed_strike;
    };

    // check validity of option parameters
    void check() const
    {
        assert(dates.size() > 0);
        assert(dates.size() == spot.size());
        for(size_t i=0; i<dates.size()-1; i++) {
            assert(dates[i]<dates[i+1]);
            if(dates[i]<=t) {
                assert(spot[i]>0.0);
            }
        }
        assert(dates.back()<=T);
    }

    // payoff profile of the option at maturity T
    double payoff(double s, double a) const
    {
        check();
        double  intrinsic=0.0;
        double  pc = (putcall==type::call) ? 1.0 : -1.0;
        switch(fixed_float) {
        case type::fixed_strike:
            intrinsic=std::max(pc*(a-strike),0.0);
            break;
        case type::floating_strike:
            intrinsic=std::max(pc*(s-strike*a),0.0);
            break;
        default:
            // option payoff not implemented
            assert(false);
        }
        return intrinsic;
    }
    void current_avg(double& avg, int& pastobs, int& obs) const
    {
        assert(dates.size()==spot.size());
        avg=0.0;
        size_t i=0;
        while(i<dates.size() && dates.at(i)<=t) {
            avg+=spot[i];
            i++;
        }
        if(i>0)
            avg /= (double) i;
        pastobs=i;
        obs=(int) dates.size();
    }


};


} // namespace fin

#endif // MF_FDM_OPTION_HPP
