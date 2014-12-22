#pragma once

#include "COMMON_DEFINITIONS.h"

template<class TT>
class ARRAY_3D
{
public: // Essential Data
	union
	{
		struct{int i_start, j_start, k_start, i_end, j_end, k_end;};
		struct{int ix_start[3], ix_end[3];};
	};

public: // Data Array
	TT* values;

public: // Speedup Variables
	int i_res, j_res, k_res;
	int ij_res, ijk_res;

public: // Constructors and Destructor
	ARRAY_3D(void)
		: values(0)
	{}

	ARRAY_3D(const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const bool& initialize = false)
		: values(0)
	{
		Initialize(i_start_input, j_start_input, k_start_input, i_res_input, j_res_input, k_res_input, initialize);
	}

	ARRAY_3D(const VI& start_ix_input, const VI& res_ix_input, const bool& initialize = false)
		: values(0)
	{
		Initialize(start_ix_input, res_ix_input, initialize);
	}

	~ARRAY_3D(void)
	{
		if(values != 0) delete [] values;
	}

public: // Initialization Function
	void Initialize(const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const bool& initialize = false)
	{
		if(values != 0) delete [] values;

		i_start = i_start_input;
		j_start = j_start_input;
		k_start = k_start_input;

		i_res = i_res_input;
		j_res = j_res_input;
		k_res = k_res_input;

		i_end = i_start + i_res - 1;
		j_end = j_start + j_res - 1;
		k_end = k_start + k_res - 1;

		ij_res = i_res*j_res;
		ijk_res = ij_res*k_res;

		assert(i_res > 0 && j_res > 0 && k_res > 0);
		values = new TT [ijk_res];

		if(initialize == true) AssignAllValues(TT());
	}

	void Initialize(const VI& start_ix_input, const VI& res_ix_input, const bool& initialize = false)
	{
		Initialize(start_ix_input.i, start_ix_input.j, start_ix_input.k, res_ix_input.i, res_ix_input.j, res_ix_input.k, initialize);
	}

public: // Operator Overloading
	TT& operator[](const int& ix) const
	{
		assert(ix >= 0 && ix <= ijk_res);

		return values[ix];
	}
	
	TT& operator()(const int& ix) const
	{
		assert(ix >= 0 && ix <= ijk_res);

		return values[ix];
	}

	TT& operator()(const int& i, const int& j, const int& k) const
	{
		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);
		assert(k >= k_start && k <= k_end);

		// TODO: Check performance 1. Use Index Function Call, 2. Use Pointer Operation for Indexing
		// return values[Index1D(i,j,k)];
    	return *(values + (i - i_start) + (j - j_start)*i_res + (k - k_start)*ij_res);
	}

	TT& operator()(const VI& ix) const
	{
		return (*this)(ix.i, ix.j, ix.k);
	}

	void operator*=(const T& constant)
	{
		for(int i = 0; i < ijk_res; i++)
		{
			values[i] *= constant;
		}
	}

	void operator +=(const T& constant)
	{
		for(int i = 0; i < ijk_res; i++)
		{
			values[i] += constant;
		}
	}

	void operator -=(const T& costant)
	{
		for(int i = 0; i < ijk_res; i++)
		{
			values[i] -= constant;
		}
	}

public: // Indexing Functions
	inline const int Index1D(const int& i, const int& j, const int& k) const
	{
		// Always use this when you program something!! It prevent from fatal errors
		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);
		assert(k >= k_start && k <= k_end);

		return (i - i_start) + (j - j_start)*i_res + (k - k_start)*ij_res;
	}

	inline const int Index1D(const VI& index) const
	{
		assert(index.i >= i_start && index.i <= i_end);
		assert(index.j >= j_start && index.j <= j_end);
		assert(index.k >= k_start && index.k <= k_end);

		return (index.i - i_start) + (index.j - j_start)*i_res + (index.k - k_start)*ij_res;
	}

	const VI Get3DIndex(const int& index_1d) const
	{
		const int k = index_1d / ij_res; 
		const int j = (index_1d - k*ij_res) / j_res;
		const int i = index_1d - k*ij_res - j*i_res;
		
		return VI(i, j, k);
	}
	
	inline const VI ClampedIndex(const int& i, const int& j, const int& k) const
	{
		return VI(CLAMP(i,i_start,i_end), CLAMP(j,j_start,j_end), CLAMP(k,k_start,k_end));
	}

	inline const VI ClampedIndex(const VI& ix) const
	{
		return ClampedIndex(ix.i, ix.j, ix.k);
	}

	inline bool Inside(const VI& ix) const
	{
		if(ix.i < i_start) return false;
		else if(ix.i > i_end) return false;
		else if(ix.j < j_start) return false;
		else if(ix.j > j_end) return false;
		else if(ix.k < k_start) return false;
		else if(ix.k > k_end) return false;
		else return true;
	}

	inline bool Inside(const int& i, const int& j, const int& k) const
	{
		if(i < i_start) return false;
		else if(i > i_end) return false;
		else if(j < j_start) return false;
		else if(j > j_end) return false;
		else if(k < k_start) return false;
		else if(k > k_end) return false;
		else return true;
	}

public: // Member Functions
	TT& DeviatedX(const int& base, const int& i_deviation) const
	{
		assert((base + i_deviation) >= 0 && (base + i_deviation) <= ijk_res);

		return values[base + i_deviation];
	}

	TT& DeviatedY(const int& base, const int& j_deviation) const
	{
		assert((base + j_deviation*i_res) >= 0 && (base + j_deviation*j_res) <= ijk_res);

		return values[base + j_deviation*i_res];
	}

	TT& DeviatedZ(const int& base, const int& k_deviation) const
	{
		asset((base + k_deviation*ij_res) >= 0 && (base + k_deviation*ij_res) <= ijk_res);

		return values[base + k_deviation*ij_res];
	}
		
	void AssignAllValues(const TT& constant)
	{
		for(int i = 0; i < ijk_res; i++)
		{
			values[i] = constant;
		}
	}

	void AssignRegionalValues(const TT& constant, const int& i_start, const int& j_start, const int& k_start, const int& i_end, const int& j_end, const int& k_end)
	{
		for(int i = i_start; i <= i_end; i++)
		{
			for(int j = j_start; j <= j_end; j++)
			{
				for(int k = k_start; k <= k_end; k++)
				{
					values[Index1D(i,j,k)] = constant;
				}
			}
		}
	}
	
	inline void Assign(const int& i, const int& j, const int& k, const TT& value)
	{
		*(values + (i - i_start) + (j - j_start)*i_res + (k - k_start)*ij_res) = values;
	}
		
// These will be added later
public: // Read and Write Functions
	void Read()
	{}

	void Write()
	{}
};








