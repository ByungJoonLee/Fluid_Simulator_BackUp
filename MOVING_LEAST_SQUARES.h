#ifndef __MOVING_LEAST_SQUARES__
#define __MOVING_LEAST_SQUARES__

#include "MATRIX_MXN.h"
#include "MATRIX_NXN.h"
#include "COMMON_DEFINITIONS.h"
#include "VECTOR_ND.h"
#include "VECTOR_3D.h"
#include "ARRAY.h"
#include <vector>

#define MLS_VALUE						0
#define MLS_VALUE_AND_NORMAL			1
#define MLS_DIFFUSIVE_NORMAL_COMPONENT	2
#define MLS_DIFFUSIVE_NORMAL			3

class MOVING_LEAST_SQUARES
{
public: // Generic Data
	const static int dimension = 3;
	bool delete_constraints_when_destroyed;

public:
	class MLS_CONSTRAINT
	{
	public: // Essential Data
		int			type;
		T			value, weight, divergence;
		VT			position, normal;
		ARRAY<VT>	vector_gradient;
		T			epsilon;

	public: // Constructors and Destructor
		MLS_CONSTRAINT(const VT& position_input, const T& value_input)
		{
			position = position_input;
			value = value_input;
			type = 0;
			weight = 1;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const VT& position_input, const T& value_input, const int& type_input)
		{
			position = position_input;
			value = value_input;
			type = type_input;
			weight = 1;
			epsilon = -1;
		}
		
		MLS_CONSTRAINT(const VT& position_input, const T& value_input, const VT& normal_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = 1;
			weight = 1;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = 1;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input, const T& weight_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = weight_input;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input, const T& weight_input, const T& divergence_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = weight_input;
			divergence = divergence_input;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input, const T& weight_input, const T& divergence_input, const ARRAY<VT>& vector_gradient_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = weight_input;
			divergence = divergence_input;
			vector_gradient = vector_gradient_input;
			epsilon = -1;
		}

		~MLS_CONSTRAINT(void)
		{}

	public: // Member Functions
		void Set(const VT& position_input, const T& value_input)
		{
			position = position_input;
			value = value_input;
			type = 0;
			weight = 1;
		}

		void Set(const VT& position_input, const T& value_input, const int& type_input)
		{
			position = position_input;
			value = value_input;
			type = type_input;
			weight = 1;
		}

		void Set(const VT& position_input, const T& value_input, const VT& normal_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = 1;
			weight = 1;
		}
		
		void Set(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = 1;
		}

		void Set(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input, const T& weight_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = weight_input;
		}

		void Set(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input, const T& weight_input, const T& divergence_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = weight_input;
			divergence = divergence_input;
		}

		void Set(const VT& position_input, const T& value_input, const VT& normal_input, const int& type_input, const T& weight_input, const T& divergence_input, const ARRAY<VT>& vector_gradient_input)
		{
			position = position_input;
			value = value_input;
			normal = normal_input;
			type = type_input;
			weight = weight_input;
			divergence = divergence_input;
			vector_gradient = vector_gradient_input;
		}
	};

public: // Essential Data
	T epsilon, epsilon_squared, one_over_epsilon_squared;
	int size_of_basis, degree_of_basis;

public: // Constraints
	vector<MLS_CONSTRAINT*> constraints;

public: // Constructors and Destructor
	MOVING_LEAST_SQUARES(void)
	{}

	MOVING_LEAST_SQUARES(const int& degree_of_basis_input = 0, const T& epsilon_input = (T)0.001, const int& num_of_constraints_input = 0)
	{
		Initialize(degree_of_basis_input, epsilon_input, num_of_constraints_input);
	}

	~MOVING_LEAST_SQUARES()
	{
		if(delete_constraints_when_destroyed == true)
		{
			const int number_of_constraints = (int)constraints.size();
			for(int i = 0; i < number_of_constraints; i++)
				delete constraints[i];
		}
	}

public: // Initialization Functions
	void Initialize(const int& degree_of_basis_input, const T& epsilon_input, const int& number_of_constraints_input)
	{
		delete_constraints_when_destroyed = true;

		SetDegreeOfBasis(degree_of_basis_input);
		SetEpsilon(epsilon_input);

		constraints.reserve(number_of_constraints_input);
		for(int i = 0; i < number_of_constraints_input; i++) constraints.push_back(new MLS_CONSTRAINT(VT(), (T)0));
	}

public: // Member Functions
	inline void Reset()
	{
		int number_of_constraints = (int)constraints.size();
		for(int i = 0; i < number_of_constraints; i++)
			delete constraints[i];

		constraints.clear();
	}

	inline void AddConstraint(const VT& position, const T& value)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value));
	}

	inline void AddConstraint(const VT& position, const T& value, const int& type)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, type));
	}

	inline void AddConstraint(const VT& position, const T& value, const VT& normal)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal));
	}

	inline void AddConstraint(const VT& position, const T& value, const VT& normal, const int& type)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal, type));
	}

	inline void AddConstraint(const VT& position, const T& value, const VT& normal, const int& type, const T& weight)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal, type, weight));
	}

	inline void AddConstraint(const VT& position, const T& value, const VT& normal, const int& type, const T& weight, const T& divergence)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal, type, weight, divergence));
	}

	inline void AddConstraint(const VT& position, const T& value, const VT& normal, const int& type, const T& weight, const T& divergence, const ARRAY<VT>& vector_gradient)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal, type, weight, divergence, vector_gradient));
	}

	inline void AddConstraint(MLS_CONSTRAINT* mls_constraint)
	{
		constraints.push_back(mls_constraint);
	}

	void SetDegreeOfBasis(const int& degree_of_basis_input)
	{
		degree_of_basis = degree_of_basis_input;

		if(degree_of_basis == 0) size_of_basis = 1;
		else if(degree_of_basis == 1)
			size_of_basis = 1 + dimension;
		else if(degree_of_basis == 2)
		{
			if(dimension == 2) size_of_basis = 6;
			else if(dimension == 3) size_of_basis = 10;
		}
		else if(degree_of_basis == 3)
		{
			if(dimension == 2) size_of_basis = 10;
			else if(dimension == 3) size_of_basis = 20;
		}
		else
		{
			cout << "MLS degree of basis is not defined with" << degree_of_basis << "degree" << endl;
			exit(1);
		}
	}
	
	void SetEpsilon(const T& epsilon_input)
	{
		epsilon = epsilon_input;
		epsilon_squared = epsilon*epsilon;
		one_over_epsilon_squared = (T)1/epsilon_squared;
	}

	VECTOR_ND<T> GetBasis(const VT& position)
	{
		VECTOR_ND<T> basis;
		basis.Initialize(size_of_basis, true);
		if(degree_of_basis == 0) basis[0] = (T)1;
		else if(degree_of_basis = 1)
		{
			basis[0] = (T)1;
			for(int i = 1; i < size_of_basis; i++) basis[i] = position.values[i-1];
		}
		else if(degree_of_basis == 2)
		{
			if(dimension == 2)
			{
				T x = position.values[0], y = position.values[0];
				basis[0] = 1;
				basis[1] = x;
				basis[2] = y;
				basis[3] = x*x;
				basis[4] = x*y;
				basis[5] = y*y;
			}
			else if(dimension == 3)
			{
				T x = position.values[0], y = position.values[1], z = position.values[2];
				basis[0] = 1;
				basis[1] = x;
				basis[2] = y;
				basis[3] = z;
				basis[4] = x*x;
				basis[5] = y*y;
				basis[6] = z*z;
				basis[7] = x*y;
				basis[8] = y*z;
				basis[9] = x*z;
			}
		}
		else if(degree_of_basis == 3)
		{
			if(dimension == 2)
			{
				T x = position.values[0], y = position.values[1];
				basis[0] = 1;
				basis[1] = x;
				basis[2] = y;
				basis[3] = x*x;
				basis[4] = y*y;
				basis[5] = x*y;
				basis[6] = x*y*y;
				basis[7] = x*x*y;
				basis[8] = x*x*x;
				basis[9] = y*y*y;
			}
			else if(dimension == 3)
			{
				T x = position.values[0], y = position.values[1], z = position.values[2];
				basis[0] = 1;
				basis[1] = x;
				basis[2] = y;
				basis[3] = z;
				basis[4] = x*x;
				basis[5] = y*y;
				basis[6] = z*z;
				basis[7] = x*y;
				basis[8] = y*z;
				basis[9] = x*z;
				basis[10] = x*x*y;
				basis[11] = x*x*z;
				basis[12] = x*y*y;
				basis[13] = y*y*z;
				basis[14] = x*z*z;
				basis[15] = y*z*z;
				basis[16] = x*y*z;
				basis[17] = x*x*x;
				basis[18] = y*y*y;
				basis[19] = z*z*z;
			}	
		}

		return basis;
	}

	VECTOR_ND<T> GetBasisDerivative(const VT& position, const int& d)
	{
		VECTOR_ND<T> basis_derivative;
		basis_derivative.Initialize(size_of_basis, true);

		if(degree_of_basis == 0)
			basis_derivative[0] = (T)0;
		else if(degree_of_basis == 1)
		{
			basis_derivative[0] = (T)0;
			if(size_of_basis > 1)
				basis_derivative[d+1] = (T)1;
		}
		else if(degree_of_basis == 2)
		{
			if(dimension == 2)
			{
				T x = position.values[0],y = position.values[1];

				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 1;
					basis_derivative[2] = 0;
					basis_derivative[3] = 2 * x;
					basis_derivative[4] = 0;
					basis_derivative[5] = y;
				}
				else if(d == 1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 1;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2 * y;
					basis_derivative[5] = x;
				}
			}
			else if(dimension == 3)
			{
				T x=position.values[0], y = position.values[1], z = position.values[2];
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 1;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2 * x;
					basis_derivative[5] = 0;
					basis_derivative[6] = 0;
					basis_derivative[7] = y;
					basis_derivative[8] = 0;
					basis_derivative[9] = z;
				}
				else if(d == 1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 1;
					basis_derivative[3] = 0;
					basis_derivative[4] = 0;
					basis_derivative[5] = 2 * y;
					basis_derivative[6] = 0;
					basis_derivative[7] = x;
					basis_derivative[8] = z;
					basis_derivative[9] = 0;
				}
				else if(d == 2)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 1;
					basis_derivative[4] = 0;
					basis_derivative[5] = 0;
					basis_derivative[6] = 2 * z;
					basis_derivative[7] = 0;
					basis_derivative[8] = y;
					basis_derivative[9] = x;
				}
			}
		}
		else if(degree_of_basis == 3)
		{
			if(dimension == 2)
			{
				T x=position.values[0], y = position.values[1];
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 1;
					basis_derivative[2] = 0;
					basis_derivative[3] = 2 * x;
					basis_derivative[4] = 0;
					basis_derivative[5] = y;
					basis_derivative[6] = y*y;
					basis_derivative[7] = 2*x*y;
					basis_derivative[8] = 3*x*x;
					basis_derivative[9] = 0;
				}
				else if(d==1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 1;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2*y;
					basis_derivative[5] = x;
					basis_derivative[6] = 2*x*y;
					basis_derivative[7] = x*x;
					basis_derivative[8] = 0;
					basis_derivative[9] = 3 * y * y;
				}
			}
			else if(dimension == 3)
			{
				T x = position.values[0], y = position.values[1], z = position.values[2];
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 1;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2 * x;
					basis_derivative[5] = 0;
					basis_derivative[6] = 0;
					basis_derivative[7] = y;
					basis_derivative[8] = 0;
					basis_derivative[9] = x;
					basis_derivative[10] = 2 * x * y;
					basis_derivative[11] = 2 * x * z;
					basis_derivative[12] = y * y;
					basis_derivative[13] = 0;
					basis_derivative[14] = z * z;
					basis_derivative[15] = 0;
					basis_derivative[16] = y * z;
					basis_derivative[17]= 3 * x * x;
					basis_derivative[18]= 0;
					basis_derivative[19]= 0;
				}
				else if(d==1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 1;
					basis_derivative[3] = 0;
					basis_derivative[4] = 0;
					basis_derivative[5] = 2 * y;
					basis_derivative[6] = 0;
					basis_derivative[7] = x;
					basis_derivative[8] = y;
					basis_derivative[9] = 0;
					basis_derivative[10] = x * x;
					basis_derivative[11] = 0;
					basis_derivative[12] = 2 * x * y;
					basis_derivative[13] = 2 * y * z;
					basis_derivative[14] = 0;
					basis_derivative[15] = z * z;
					basis_derivative[16] = x * z;
					basis_derivative[17] = 0;
					basis_derivative[18] = 3 * y * y;
					basis_derivative[19] = 0;
				}
				else if(d==2)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 1;
					basis_derivative[4] = 0;
					basis_derivative[5] = 0;
					basis_derivative[6] = 2 * z;
					basis_derivative[7] = 0;
					basis_derivative[8] = y;
					basis_derivative[9] = x;
					basis_derivative[10] = 0;
					basis_derivative[11] = x * x;
					basis_derivative[12] = 0;
					basis_derivative[13] = y * y;
					basis_derivative[14] = 2 * x * z;
					basis_derivative[15] = 2 * y * z;
					basis_derivative[16] = x * y;
					basis_derivative[17] = 0;
					basis_derivative[18] = 0;
					basis_derivative[19] = 3 * z * z;
				}
			}
		}

		return basis_derivative;
	}

	VECTOR_ND<T> GetBasisSecondDerivative(const VT& position, const int& d)
	{
		VECTOR_ND<T> basis_derivative(size_of_basis);
		
		if(degree_of_basis == 0)
			basis_derivative[0] = (T)0;
		else if(degree_of_basis == 1)
		{
			basis_derivative[0] = (T)0;
			if(size_of_basis > 1)
				basis_derivative[d+1] = (T)0;
		}
		else if(degree_of_basis == 2)
		{
			if(dimension == 2)
			{
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 2;
					basis_derivative[4] = 0;
					basis_derivative[5] = 0;
				}
				else if(d == 1)
				{					
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2;
					basis_derivative[5] = 0;
				}
			}
			else if(dimension == 3)
			{
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2;
					basis_derivative[5] = 0;
					basis_derivative[6] = 0;
					basis_derivative[7] = 0;
					basis_derivative[8] = 0;
					basis_derivative[9] = 0;
				}
				else if(d == 1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 0;
					basis_derivative[5] = 2;
					basis_derivative[6] = 0;
					basis_derivative[7] = 0;
					basis_derivative[8] = 0;
					basis_derivative[9] = 0;
				}
				else if(d==2)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 0;
					basis_derivative[5] = 0;
					basis_derivative[6] = 2;
					basis_derivative[7] = 0;
					basis_derivative[8] = 0;
					basis_derivative[9] = 0;
				}
			}
		}
		else if(degree_of_basis == 3)
		{
			if(dimension == 2)
			{
				T x = position.values[0], y = position.values[1];
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 2;
					basis_derivative[4] = 0;
					basis_derivative[5] = 0;
					basis_derivative[6] = 0;
					basis_derivative[7] = 2 * y;
					basis_derivative[8] = 6 * x;
					basis_derivative[9] = 0;
				}
				else if(d == 1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2;
					basis_derivative[5] = 0;
					basis_derivative[6] = 2 * x;
					basis_derivative[7] = 0;
					basis_derivative[8] = 0;
					basis_derivative[9] = 6 * y;
				}
			}
			else if(dimension==3)
			{
				T x = position.values[0], y = position.values[1], z = position.values[2];
				if(d == 0)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 2;
					basis_derivative[5] = 0;
					basis_derivative[6] = 0;
					basis_derivative[7] = 0;
					basis_derivative[8] = 0;
					basis_derivative[9] = 1;
					basis_derivative[10] = 2 * y;
					basis_derivative[11] = 2 * z;
					basis_derivative[12] = 0;
					basis_derivative[13] = 0;
					basis_derivative[14] =0;
					basis_derivative[15] = 0;
					basis_derivative[16] = 0;
					basis_derivative[17] = 6 * x;
					basis_derivative[18] = 0;
					basis_derivative[19] = 0;
				}
				else if(d == 1)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0; 
					basis_derivative[4] = 0;
					basis_derivative[5] = 2;
					basis_derivative[6] = 0;
					basis_derivative[7] = 0;
					basis_derivative[8] = 1;
					basis_derivative[9] = 0;
					basis_derivative[10] = 0;
					basis_derivative[11] = 0;
					basis_derivative[12] = 2 * x;
					basis_derivative[13] = 2 * z;
					basis_derivative[14] = 0;
					basis_derivative[15] = 0;
					basis_derivative[16] = 0;
					basis_derivative[17] = 0;
					basis_derivative[18] = 6 * y;
					basis_derivative[19] = 0;
				}
				else if(d == 2)
				{
					basis_derivative[0] = 0;
					basis_derivative[1] = 0;
					basis_derivative[2] = 0;
					basis_derivative[3] = 0;
					basis_derivative[4] = 0;
					basis_derivative[5] = 0;
					basis_derivative[6] = 2;
					basis_derivative[7] = 0;
					basis_derivative[8] = 0;
					basis_derivative[9] = 0;
					basis_derivative[10] = 0;
					basis_derivative[11] = 0;
					basis_derivative[12] = 0;
					basis_derivative[13] = 0;
					basis_derivative[14] = 2 * x;
					basis_derivative[15] = 2 * y;
					basis_derivative[16] = 0;
					basis_derivative[17] = 0; 
					basis_derivative[18] = 0;
					basis_derivative[19] = 6 * z;
				}
			}
		}
		return basis_derivative;
	}

	// Get Coefficients of interpolation polynomial
	T GetScalar(const VT& x)
	{
		const int number_of_constraints = (int)constraints.size();

		int number_of_rows = 0;
		for(int i = 0; i < number_of_constraints; i++)
		{
			switch(constraints[i]->type)
			{
			case MLS_VALUE:
				number_of_rows++;
				break;
			case MLS_VALUE_AND_NORMAL:
				number_of_rows++;
				break;
			case MLS_DIFFUSIVE_NORMAL:
				number_of_rows += dimension;
				break;
			}
		}

		MATRIX_MXN<T> B(number_of_rows, size_of_basis);
		VECTOR_ND<T> W(number_of_rows), WWphi(number_of_rows);

		int row = 0;
		for(int i = 0; i < number_of_constraints; i++)
		{
			MLS_CONSTRAINT* ctr = constraints[i];
			assert(ctr->epsilon > 0);

			VT deviation = ctr->position - x;
			// General Weight Function in MLS(w(|r|) = 1/(|r|^2 + epsilon^2)
			T w = (T)1/(deviation.SqrMagnitude() + ctr->epsilon*ctr->epsilon)*ctr->weight;
			
			if(ctr->type == MLS_VALUE)
			{
				W.values[row] = w;
				B.Set_Row(row, GetBasis(deviation));

				WWphi.values[row] = w*w*ctr->value;

				row++;
			}
			else if(ctr->type == MLS_VALUE_AND_NORMAL)
			{
				W.values[row] = w;
				B.Set_Row(row, GetBasis(deviation));

				WWphi.values[row] = w*w*(ctr->value - DotProduct(deviation, ctr->normal));

				row++;
			}
			else if(ctr->type == MLS_DIFFUSIVE_NORMAL)
			{
				for(int d = 0; d < dimension; d++)
				{
					W.values[row] = w;
					VECTOR_ND<T> value = GetBasisDerivative(deviation, d);

					B.Set_Row(row, value);

					WWphi.values[row] = w*w*ctr->normal.values[d];
					row++;
				}
			}
		}

		// Solve the normal equation for calculating the coefficients
		VECTOR_ND<T> c = (B.Weighted_Normal_Equations_Matrix(W)).Cholesky_Solve(B.Transposed_Multiply(WWphi));

		return c.values[0];
	}
};
#endif
