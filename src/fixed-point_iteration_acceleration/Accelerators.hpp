//! @file Accelerators.hpp
//! @brief Header file for the (polymorphic family of) iterators

#ifndef SRC_ACCELERATORS_HPP_
	#define SRC_ACCELERATORS_HPP_
	
	#include <deque>
	#include <cassert>
	#include <algorithm>
	#include "FixedPointTraits.hpp"
	#include "RotatingMatrix.hpp"
	
	#ifndef Traits
		//! All the classes in this project derive from Traits to have uniform types.
		#define Traits SparseTraits
		// #define Traits DenseTraits // choose this or the previous one (you can also do it in the Makefile)
	#endif
	namespace FixedPoint
	{
		//! @brief A function object producing iterations of a sequence.
		//!
		//! This class is the parent class for the subsequent classes: what it does
		//! is providing, through Iterator::operator()(const std::deque < Vector > & past),	
      //! the next value of a (hopefully converging) sequence, making use of an
		//! iteration function #phi and previous iterations data (only the last value in the base case).
		class Iterator : public Traits
		{
			public:
			
			//! @brief  The constructor
			//!
			//! @tparam IterationFun   An IterationFunction type
			//! @param  IF             An IterationFunction to forward to the class
			//! @param  dim            A std::size_t representing the dimension of the space in which the IterationFunction operates
			//! @sa dimension DenseTraits SparseTraits
			template <class IterationFun>
			Iterator(IterationFun && IF , std::size_t dim): phi(std::forward<IterationFun>(IF)), dimension (dim) {}
			
			//! @brief Function call operator
			//!
			//! @param  past  The sequence of past iterates used to compute a new value
			//! @return       The new computed value
			virtual inline Vector operator()(const std::deque < Vector > & past); // definition is at the end of this file
			
			//! Resetting the Iterator. The meaning will be clear in the derived classes.
			virtual void reset(){}
			
			//! SetUp method. For AndersonAccelerator class.
			virtual void SetUp (const std::deque<Vector> &){}
			
			//! Destructor is virtual for the classes that will inherit from here
			virtual ~Iterator(){};
			
			//! Get #phi
			IterationFunction & getIterationFunction(){ return phi;}
			//! Get #phi (const version)
			IterationFunction getIterationFunction() const { return phi;}
			//! Get #dimension
			std::size_t & getDimension() { return dimension ;}
			//! Get #dimension (const version)
			std::size_t getDimension() const { return dimension ;}
			
			
			protected:
			//! A function computing a new vector from an older one
			IterationFunction phi;
			//! Fixed-point problems go from \f$\mathbb{R}^{dimension}\f$ to \f$\mathbb{R}^{dimension}\f$
			std::size_t dimension;
		};
		
		//! @brief Alternate secant method
		//!
		//! This is a two-level acceleration method equivalent to Anderson with memory two.
		//! \f[ x_{n+1}=\phi(x_n)-\frac{(\Delta x_n -\Delta x_{n-1})\cdot\Delta x_n}{||\Delta x_n -\Delta x_{n-1}||^2}(\phi(x_n)-\phi(x_{n-1}) \f]
		//!
		//! Reference: V. Eyert, A Comparative Study on Methods for Convergence Acceleration of Iterative Vector Sequences, Journal of Computational Physics 124 (1996) 271–285.
		//! H. Fang, Y. Saad, Two classes of multisecant methods for nonlinear acceleration, Numerical Linear Algebra with Applications 16 (2009) 197–221.		
		class ASecantAccelerator : public Iterator
		{
			public:
			//! @brief Constructor
			//!
			//! Delegates to the parent class constructor and allocates memory for the vectors.
			//! @sa Iterator::Iterator
			template <class IterationFun>
			ASecantAccelerator(IterationFun&& IF, std::size_t dim):
			Iterator(std::forward<IterationFun>(IF), dim), deltaXOld(dim), phiOld(dim), firstTime(true) {}
			
			//! Update formula is in the description of class ASecantAccelerator
			Vector operator()(const std::deque < Vector > &) override;
			
			//! Resets the #firstTime bool
			void reset() override {firstTime=true;}
			
			private:
			
			//! Being a two-level algorithm, I need keep in memory \f$\Delta x_{n-1}\f$
			Vector deltaXOld;
			//! Being a two-level algorithm, I need keep in memory \f$\phi(x_{n-1})\f$
			Vector phiOld;
			//! Indicates if the call operator has never been called or not
			bool firstTime;
		};
		
		//! @brief Anderson Accelerator
		//!
		//! \f[x_{k+1}=x_{k} + \beta f_k -(\mathscr{X}_k + \beta \mathscr{F}_k)(\mathscr{F}_k^T \mathscr{F}_k)^{-1}\mathscr{F}_k^T f_k \f]
		//! where
		//!
		//! \f$\mathscr{X}_k=[\Delta x_{k-m}...\Delta x_{k-1}]\f$ with \f$\Delta x_{i}=x_{i+1}-x_i\f$
		//!
		//! \f$\mathscr{F}_k=[\Delta f_{k-m}...\Delta f_{k-1}]\f$ with \f$\Delta f_{i}=f_{i+1}-f_i\f$ and \f$f_i\f$ is \f$\phi(x_i)-x_{i}\f$
		//!
		//! Reference: Anderson, Donald, Iterative Procedures for Nonlinear Integral Equations J. ACM 12 (1965) 547-560, doi 10.1145/321296.321305
		class AndersonAccelerator : public Iterator
		{
			public:
			//! @brief Constructor
			//!
			//! Delegates to the parent class constructor and sets #mixingParameter and #memory
			//! @sa Iterator::Iterator
			template <class IterationFun>
			AndersonAccelerator(IterationFun&& IF, std::size_t dim, double p =1, std::size_t m = 10): 
			Iterator(std::forward<IterationFun>(IF), dim ), mixingParameter(p), memory(std::min (m, dim)) {}
			
			//! @brief Update formula is in the description of class AndersonAccelerator
			//!
			//! Very important: the Matrices #X and #F must be in the correct state, otherwise call first AndersonAccelerator::SetUp
			Vector operator()(const std::deque < Vector > &) override;
			
			//! Resets #X and #F
			void reset() override 
			{
				X.reset();
				F.reset();
			}
			
			//! @brief For advanced use
			//!
			//! In the case AndersonAccelerator is called with \p past that was not built up by AndersonAccelerator itself,
			//! this method sets matrices #X and #F and vector #fOld in the proper state
			void SetUp (const std::deque < Vector > & past) override
			{
				reset();
				if ( past.size() < 2 or memory < 2 ) return;
				if ( past.size() == 2 or memory == 2) {
					fOld = phi(past[past.size()-2]) - past[past.size()-2];
				}
				else 
				{
					int m_k = std::min( memory, past.size() )- 1u;
					fOld = phi(past[past.size()-2]) - past[past.size()-2];
					for (int j = 0; j < m_k -1u ; ++j )
					{
						X.push_back( past [ past.size()-m_k+j ]- past [ past.size()-m_k+j-1 ] );		
					}
					
					for (int j = 0; j < m_k - 1u; ++j )
					{
						F.push_back
						(
							phi ( past [ past.size()-m_k+j ] ) -
							phi ( past [ past.size()-m_k+j-1 ] ) -
							X.col(j)
						);
					}		
				}
			}
			private:
			//! A sort of relaxation parameter.
			double mixingParameter;
			
			//! How many past iterates are considered at each step.
			std::size_t memory;
			//!previous residual vector
			Vector fOld;
			
			
			//! See update formula in AndersonAccelerator
			apsc::RotatingMatrixX <double , apsc::InsertStrategy::NewestReplacesOldest > X {dimension, memory -1};
			//! See update formula in AndersonAccelerator
			apsc::RotatingMatrixX <double , apsc::InsertStrategy::NewestReplacesOldest > F { dimension, memory -1};
		};
		
		//! The simplest iterator: evaluates #phi in the last value of \p past 
		Traits::Vector Iterator::operator()(const std::deque < Vector > & past)
		{
			assert (!past.empty());
			return phi( past.back() );
		}
		
	}
#endif /* SRC_ACCELERATORS_HPP_ */					