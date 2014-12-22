#pragma once

#include "COMMON_DEFINITIONS.h"
#include "MULTITHREADING.h"

template<class TT>
class VECTOR_ND
{
public: // Essential Data
	int		num_dimension;
	TT*		values;
	
public: // Constructors and Destructor
	VECTOR_ND(void)
		: num_dimension(0), values(0)
	{}

	VECTOR_ND(const int& num_dimension_input)
	{
		values = 0;
		
		Initialize(num_dimension_input);
	}

	VECTOR_ND(const VECTOR_ND<TT>& vector_input)
	{
		values = 0;
		
		Initialize(vector_input.num_dimension, false);
		
		for (int i = 0; i < num_dimension; i++)
		{
			values[i] = vector_input[i];
		}
	}

	~VECTOR_ND(void)
	{
		if (values != 0) 
		{
			delete [] values;
		}
		
		num_dimension = 0;
	}

public: // Initialization Functions
	void Initialize(const int& num_dimension_input, const bool initialize = false)
	{
		num_dimension = num_dimension_input;

		DELETE_POINTER(values);
		
		if (num_dimension > 0)
		{
			values = new TT [num_dimension];
						
			if(initialize == true)
			{
				for(int i = 0; i < num_dimension; i++) 
				{
					values[i] = TT();
				}
			}
		}
	}

public: // Operator Overloading
	inline void operator = (const VECTOR_ND<TT>& from)
	{
		num_dimension = from.num_dimension;

		DELETE_POINTER(values);
		
		values = new TT [num_dimension];
		
		for(int i = 0; i < num_dimension; i++) 
		{
			values[i] = from[i];
		}
	}

	inline TT& operator[](const int& i) const
	{
		return values[i];
	}

	inline VECTOR_ND<TT> operator+(const VECTOR_ND<TT>& vector) const
	{
		assert(num_dimension == vector.num_dimension);

		VECTOR_ND<TT> Add_result(num_dimension);
		
		for(int i = 0; i < num_dimension; i++) Add_result.values[i] = values[i] + vector.values[i];

		return Add_result;
	}

	inline VECTOR_ND<TT> operator-(const VECTOR_ND<TT>& vector) const
	{
		assert(num_dimension == vector.num_dimension);

		VECTOR_ND<TT> Subtract_result(num_dimension);

		for(int i = 0; i < num_dimension; i++) Subtract_result[i] = values[i] - vector.values[i];

		return Subtract_result;
	}

	void operator += (const TT& s)
	{
		for(int i = 0; i < num_dimension; i++) values[i] += s;
	}

	void operator -= (const TT& s)
	{
		for(int i = 0; i < num_dimension; i++) values[i] -= s;
	}

	void operator *= (const int& s)
	{
		for(int i = 0; i < num_dimension; i++) values[i] *= s;
	}

	void operator *= (const T& s)
	{
		for(int i = 0; i < num_dimension; i++) values[i] *= s;
	}

	void operator /= (const TT& s)
	{
		TT one_over_s = (TT)1/s;
	
		for(int i = 0; i < num_dimension; i++) values[i] *= one_over_s;
	}

	void operator += (const VECTOR_ND<TT>& s)
	{
		assert(num_dimension == s.num_dimension);

		for(int i = 0; i <num_dimension; i++) values[i] += s.values[i];
	}

	void operator -= (const VECTOR_ND<TT>& s)
	{
		assert(num_dimension == s.num_dimension);

		for(int i = 0; i < num_dimension; i++) values[i] -= s.values[i];
	}

	void operator *= (const VECTOR_ND<TT>& s)
	{
		assert(num_dimension == s.num_dimension);

		for(int i = 0; i < num_dimension; i++) values[i] *= s.values[i];
	}

	void operator /= (const VECTOR_ND<TT>& s)
	{
		assert(num_dimension == s.num_dimension);

		for(int i = 0; i < num_dimension; i++) values[i] /= s.values[i];
	}

	VECTOR_ND<TT> operator*(const T& s) const
	{
		VECTOR_ND<TT> Multiply_result(num_dimension);
		
		for(int i = 0; i < num_dimension; i++) Multiply_result.values[i] = values[i]*s;

		return Multiply_result;
	}

public: // Member Functions
	T SqrMagnitude() const
	{
		T sum(0);
		
		for(int i = 0; i < num_dimension; i++)
			sum += SQUARE(values[i]);

		return sum;
	}

	T MaxAbs() const
	{
		T max_abs(0);
		
		for(int i = 0; i < num_dimension; i++)
			max_abs = max(abs(values[i]), max_abs);
		
		return max_abs;
	}
};

inline static T DotProduct(const VECTOR_ND<T>& v1, const VECTOR_ND<T>& v2)
{
	assert(v1.num_dimension == v2.num_dimension);

	T sum = 0;
	
	for(int i = 0; i < v1.num_dimension; i++)
		sum += v1.values[i]*v2.values[i];

	return sum;
}

inline static T DotProduct(const VECTOR_ND<T>& v1, const VECTOR_ND<T>& v2, const int& start_ix, const int& end_ix)
{
	assert(v1.num_dimension == v2.num_dimension);

	T sum = 0;

	for(int i = start_ix; i <= end_ix; i++)
		sum += v1.values[i]*v2.values[i];
	
	return sum;
}

inline static T DotProduct(MULTITHREADING* multithreading, const int& thread_id, const VECTOR_ND<T>& v1, const VECTOR_ND<T>& v2)
{
	assert(v1.num_dimension == v2.num_dimension);

	T *v1_values = v1.values, *v2_values = v2.values;

	const int start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

	T dot_product = (T)0;

	for(int i = start_ix; i <= end_ix; i++) dot_product += v1_values[i]*v2_values[i];

	multithreading->SyncSum(thread_id, dot_product);

	return dot_product;
}

inline static void DotProduct(MULTITHREADING* multithreading, const int& thread_id, const VECTOR_ND<T>& v1, const VECTOR_ND<T>& v2, T& dot_result)
{
	assert(v1.num_dimension == v2.num_dimension);

	T *v1_values = v1.values, *v2_values = v2.values;

	const int start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

	dot_result = (T)0;

	for(int i = start_ix; i <= end_ix; i++) dot_result += v1_values[i]*v2_values[i];

	multithreading->SyncSum(thread_id, dot_result);
}

template<class T> inline ostream& operator << (ostream& output, const VECTOR_ND<T>& v)
{
	for(int i = 0; i < v.num_dimension; i++) output << v.values[i] << " ";
	output << endl;
	return output;
}






