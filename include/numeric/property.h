//----------------------------  property.h  ---------------------------
// 
//
//  Description:	property_type
//
//  Documents:		Property.h
//					Property.cpp 
//
//  Creation date  : 03/2009
//       modified  : 
//
//  Author: Siegel Pascal
//
//  Copyright (C) 2009 by In2Earth modelling 
//
//---------------------------- property.h  ---------------------------

#ifndef property_h
#define property_h


// system includes

#include <string>
#include <vector>

// project includes 

// local includes 

// namespace 
namespace fractures_intersect  {


/**
 *    property :	This class is for storing a property. Properties 
 *					scalars apply to  the elements of the mesh
 *
  *    @author Pascal Siegel 
 *    @date 3/2009
 */

template <typename T>  
class property
{

// variables

public:

	/**
	 *	name of the property
	 */
    std::string name_;
	
	
	/**
	 *	no data value 
	 */
	T no_data_value_;

	/**
	 *	vector of discretise property 
	 */
	
	std::vector<T> values_;
		

// functions


public : 


	/**
	 *      Default constructor 
	 */

    property<T>() {no_data_value_ = -9999;};

	/**
	 *	Copy constructor
	 * 
	 *      @param from The value to copy to this object
	 */

    property<T>(const property<T>& from)
    {
        name_ = from.name_;
        no_data_value_ = from.no_data_value_;
        values_ = from.values_;
    };

    /**
     *	constructor
     **/

     property<T>(std::string name, int size )
     {
          name_ = name;
          no_data_value_ = -9999;
          resize_to_no_data_value(size);
      };


	/**
	 *	Destructor 
	 */
        
    ~property<T>() {}; 


    /** 
     *		Get the value of the property ie element   
	 *
	 *		@param i an integer representing the desired index of element 
     * 
     *      @return  a  value _ of type T 
     */
	
    T value(int i)
    {
       return values_[i];
    };

	/** 
     *		set the value of the property ie element   
	 *
	 *		@param i an integer representing the desired index of element
	 *
	 *		@param value the value to be set  
     * 
     *      @return  a  value _ of type T 
     */
	
    void set_value(int i, T value)
    {
        values_[i]= value;
    };

	/** 
     *		set add a vclue to the property ie element   
	 *
	 *		@param i an integer representing the desired index of element
	 *
	 *		@param value the value to be set  
     * 
     *      @return  a  value _ of type T 
     */
	
    void add(int i, T value)
    {
       values_[i] += value;
    };

	/** 
     *		resize the memory for property vector   
	 *
	 *		@param i reseve size
	 */

	void resize(int i)
    {
		values_.resize(i);

        return ;
    };

	/**
     *		resize the memory for property vector
	 *
	 *		@param i reseve size
	 */

	void resize_to_no_data_value(int i)
    {
		values_.resize(i);
		init_to_no_data_value();

        return ;
    };

	/** 
     *		init to no data value
	 *
	 */

	void init_to_no_data_value()
    {
       for (int i = 0; i < values_.size(); i++) values_[i] = no_data_value_;
    };


	/**
	 *		init to zero
	 *
	 */

	void init_to_zero()
	{
		for (int i = 0; i < values_.size(); i++) values_[i] = 0;
	};

	/**
	 *		init to value
	 *
	 */

	void init_to_value(T value)
	{
		for (int i = 0; i < values_.size(); i++) values_[i] = value;
	};


	/** 
     *		reserve the memory for property vector   
	 *
	 *		@param i reseve size
	 */

	void reserve(int i)
    {
       values_.reserve(i);
    };


	/** 
     *		clear the memory for property vector   
	 *
	 */

	void clear()
    {
	   values_.clear();
       std::vector<T>().swap(values_);
    };

	/**
	 *		push back
	 *
	 */

	void push_back(T& value)
	{
		values_.push_back(value);
	};

	/**
	 *		size
	 *
	 */

	size_t size()
	{
		return values_.size();
	};

};


// inline methods


} // namespace simulator_2d


#endif  // property_h
