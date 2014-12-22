#pragma once

#include "FIELD_STRUCTURE_3D.h"
#include "DYNAMIC_ARRAY.h"

template<class TT>
class SCAN_LINE_ALGORITHM
{
public: // Essential Data
	MULTITHREADING*				multithreading;
	ARRAY<DYNAMIC_ARRAY<VI>>	stack_from_thread;

public: // Constructor and Destructor
	SCAN_LINE_ALGORITHM(void)
		: multithreading(0)
	{}

	~SCAN_LINE_ALGORITHM(void)
	{}

public: // Initialization Function
	void Initialize(MULTITHREADING* multithreading_input, const int& grid_ijk_res)
	{
		multithreading = multithreading_input;
		stack_from_thread.Initialize(multithreading->num_threads);

		int stack_size = grid_ijk_res / multithreading->num_threads;

		for (int i = 0; i < stack_from_thread.num_elements; i++)
		{
			stack_from_thread[i].Initialize(stack_size, stack_size/multithreading->num_threads);
		}
	}

public: // Member Functions
	void ClearSeedStack(const int& thread_id)
	{
		DYNAMIC_ARRAY<VI>& seed_stack = stack_from_thread[thread_id];
		seed_stack.Emptify();

		multithreading->Sync(thread_id);
	}

	void ScanFieldThreaded(FIELD_STRUCTURE_3D<TT>* field, const TT old_value, const TT new_value, const TT block_value)
	{
		multithreading->RunThreads(&SCAN_LINE_ALGORITHM::ScanField, this, field, old_value, new_value, block_value);
	}

	void ScanFieldNoBlockThreaded(FIELD_STRUCTURE_3D<TT>* field, const TT old_value, const TT new_value)
	{
		multithreading->RunThreads(&SCAN_LINE_ALGORITHM::ScanFieldNoBlock, this, field, old_value, new_value);
	}

	void ScanField(const int& thread_id, FIELD_STRUCTURE_3D<TT>* bc_field, const TT old_value, const TT new_value, const TT block_value)
	{
		FIELD_STRUCTURE_3D<TT>& field = *bc_field;

		DYNAMIC_ARRAY<VI>& seed_stack = stack_from_thread[thread_id];

		while(seed_stack.num_of_elements != 0)
			ScanLine(seed_stack.Pop(), field, thread_id, old_value, new_value, block_value);

		multithreading->Sync(thread_id);
	}

	void ScanFieldNoBlock(const int& thread_id, FIELD_STRUCTURE_3D<TT>* bc_field, const TT old_value, const TT new_value)
	{
		FIELD_STRUCTURE_3D<TT>& field = *bc_field;

		DYNAMIC_ARRAY<VI>& seed_stack = stack_from_thread[thread_id];

		while(seed_stack.num_of_elements != 0)
			ScanLine(seed_stack.Pop(), field, thread_id, old_value, new_value);
		
		multithreading->Sync(thread_id);
	}

	void ScanLine(const VI& start_idx, FIELD_STRUCTURE_3D<TT>& field, const int& current_thread_id, const TT& old_value, const TT& new_value, const TT& block_value)
	{
		GRID_STRUCTURE_3D& grid = field.grid;
		ARRAY_3D<TT>& arr = field.array_for_this;

		// Scaning process to +Y direction
		for (int j = start_idx.j; j <= grid.j_end; j++)
		{
			VI ix(start_idx.i, j, start_idx.k);
			TT& value = arr(ix);

			if(value == old_value)
			{
				value = new_value;
				AddSeed(ix, field, current_thread_id, old_value, new_value, block_value);
			}
			else
			{
				break;
			}
		}
		
		// Scaning process to -Y direction
		for (int j = start_idx.j - 1; j >= grid.j_start; j--)
		{
			VI ix(start_idx.i, j, start_idx.k);
			TT& value = arr(ix);

			if(value == old_value)
			{
				value = new_value;
				AddSeed(ix, field, current_thread_id, old_value, new_value, block_value);
			}
			else
			{
				break;
			}
		}
	}

	void AddSeed(const VI& idx, FIELD_STRUCTURE_3D<TT>& field, const int& current_thread_id, const TT& old_value, const TT& new_value, const TT& block_value)
	{
		GRID_STRUCTURE_3D& grid = field.grid;
		ARRAY_3D<TT>& arr = field.array_for_this;

		DYNAMIC_ARRAY<VI>& stack = stack_from_thread[current_thread_id];

		VI b0, b1, s;

		b0 = VI(idx.i + 1, idx.j + 1, idx.k);
		b1 = VI(idx.i + 1, idx.j - 1, idx.k);
		s  = VI(idx.i + 1, idx.j    , idx.k);
		if(arr(s) == old_value && (arr(b0) == block_value || arr(b1) == block_value))
			stack.Push(s);

		b0 = VI(idx.i - 1, idx.j + 1, idx.k);
		b1 = VI(idx.i - 1, idx.j - 1, idx.k);
		s  = VI(idx.i - 1, idx.j    , idx.k);
		if(arr(s) == old_value && (arr(b0) == block_value || arr(b1) == block_value))
			stack.Push(s);
		
		b0 = VI(idx.i, idx.j + 1, idx.k + 1);
		b1 = VI(idx.i, idx.j - 1, idx.k + 1);
		s  = VI(idx.i, idx.j    , idx.k + 1);
		if(arr(s) == old_value && (arr(b0) == block_value || arr(b1) == block_value))
			stack.Push(s);

		b0 = VI(idx.i, idx.j + 1, idx.k - 1);
		b1 = VI(idx.i, idx.j - 1, idx.k - 1);
		s  = VI(idx.i, idx.j    , idx.k - 1);
		if(arr(s) == old_value && (arr(b0) == block_value || arr(b1) == block_value))
			stack.Push(s);
	}

	void ScanLine(const VI& start_idx, FIELD_STRUCTURE_3D<TT>& field, const int& current_thread_id, const TT& old_value, const TT& new_value)
	{
		GRID_STRUCTURE_3D& grid = field.grid;
		ARRAY_3D<TT>& arr = field.array_for_this;

		// Scaning process to +Y direction
		for (int j = start_idx.j; j <= grid.j_end; j++)
		{
			VI ix(start_idx.i, j, start_idx.k);
			TT& value = arr(ix);

			if (value == old_value)
			{
				value = new_value;
				AddSeed(ix, field, current_thread_id, old_value, new_value);
			}
			else
			{
				break;
			}
		}

		// Scaning process to -Y direction
		for (int j = start_idx.j; j < grid.j_end; j++)
		{
			VI ix(start_idx.i, j, start_idx.k);
			TT& value = arr(ix);

			if (value == old_value)
			{
				value = new_value;
				AddSeed(ix, field, current_thread_id, old_value, new_value);
			}
			else
			{
				break;
			}
		}
	}

	void AddSeed(const VI& idx, FIELD_STRUCTURE_3D<TT>& field, const int& current_thread_id, const TT& old_value, const TT& new_value)
	{
		GRID_STRUCTURE_3D& grid = field.grid;
		ARRAY_3D<TT>& arr = field.array_for_this;

		DYNAMIC_ARRAY<VI>& stack = stack_from_thread[current_thread_id];

		VI b0, b1, s;

		b0 = VI(idx.i + 1, idx.j + 1, idx.k);
		b1 = VI(idx.i + 1, idx.j - 1, idx.k);
		s  = VI(idx.i + 1, idx.j    , idx.k);
		if(arr(s) == old_value && ((arr(b0) != old_value && arr(b0) != new_value) || (arr(b1) != old_value && arr(b1) != new_value)))
			stack.Push(s);

		b0 = VI(idx.i - 1, idx.j + 1, idx.k);
		b1 = VI(idx.i - 1, idx.j - 1, idx.k);
		s  = VI(idx.i - 1, idx.j    , idx.k);
		if(arr(s) == old_value && ((arr(b0) != old_value && arr(b0) != new_value) || (arr(b1) != old_value && arr(b1) != new_value)))
			stack.Push(s);

		b0 = VI(idx.i, idx.j + 1, idx.k + 1);
		b1 = VI(idx.i, idx.j - 1, idx.k + 1);
		s  = VI(idx.i, idx.j    , idx.k + 1);
		if(arr(s) == old_value && ((arr(b0) != old_value && arr(b0) != new_value) || (arr(b1) != old_value && arr(b1) != new_value)))
			stack.Push(s);

		b0 = VI(idx.i, idx.j + 1, idx.k - 1);
		b1 = VI(idx.i, idx.j - 1, idx.k - 1);
		s  = VI(idx.i, idx.j    , idx.k - 1);
		if(arr(s) == old_value && ((arr(b0) != old_value && arr(b0) != new_value) || (arr(b1) != old_value && arr(b1) != new_value)))
			stack.Push(s);
	}
};

		
