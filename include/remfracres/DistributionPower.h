#pragma once

#include <string>
#include <cmath>
#include <boost/random.hpp>

namespace remfracres {

class DistributionPower
{
  public:
	double theta_;

	double alpha_;

	DistributionPower()
	{
		theta_=1;
		alpha_= 1;
	};

	DistributionPower(double alpha, double theta)
	{
		theta_=theta;
		alpha_= alpha;
	};

    ~DistributionPower()
    {

    };
	/**
     *      Assignment operator
	 *
     *      @param from the value to assign to this object
     *
     *      @return A reference to this object
     */

    DistributionPower& operator=(const DistributionPower& from){
    	theta_ = from.theta_;
    	alpha_ = from.alpha_;
    	return *this;
    };

	/**
     *      copy function
	 *
     *      @param from the value to assign to this object
     */

	void copy(const DistributionPower& from){
    	theta_ = from.theta_;
    	alpha_ = from.alpha_;
	};

	/**
	 * @fn double GcInvPower(double)
	 * @brief
	 *
	 * @param proba
	 * @return
	 */
    double GcInvPower(boost::random::mt19937& center_gen ){
       double q;
       boost::uniform_01<boost::mt19937> center_uni(center_gen);
       double proba = center_uni();
       q = pow( proba, (1.0) / theta_) * alpha_;
       return q;
    };


};

}


