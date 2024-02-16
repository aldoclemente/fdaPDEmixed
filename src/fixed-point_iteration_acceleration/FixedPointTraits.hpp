//! @file FixedPointTraits.hpp
//! @brief Defines types used throughout the project.
//!
//! Here are defined the types that are used in the generic implementation.
//! Also some utilities related to this types are provided.
//! The project is based on Eigen, so we define a sparse and a dense implementation.
//! The user can easily change the types defined in each struct.

#ifndef _FIXEDPOINTTRAITS_HPP_
	#define _FIXEDPOINTTRAITS_HPP_
	
	#include <iostream>
	#include <Eigen/Core>
	#include <Eigen/SparseCore>
	
	//! The namespace that contains everything related to my project.
	
	namespace FixedPoint
	{
		//! @brief As the name says, uses dense matrices.
		//!
		//! Intended for small-size problems.
		struct DenseTraits
		{
			//! A dense vector.
			using Vector = Eigen::Matrix<double,Eigen::Dynamic,1>;
			//! An iteration function (in our problems are functions from \f$\mathbb{R}^n\f$ to \f$\mathbb{R}^n\f$).
			using IterationFunction = std::function <Vector (Vector const &)>;
			//! A dense matrix.
			using Matrix = Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>;
			
			//! Computes the distance between two vectors.
			static double distance (Vector const & current, Vector const & previous)
			{
				return (current-previous).norm() ;
			}
			
			//! Prints the elements of vector \p v to the output stream \p OS.
			static std::ostream & print (const Vector& v, std::ostream & OS)
			{
				for (std::size_t i = 0; i < v.size() ; ++i) OS << *(v.data()+i) << " ";
				OS << std::endl ;
				return OS;
			}
			
		};
		
		//! @brief For using sparse matrices
		//!
		//! The project is meant to be applied to large sparse matrices.
		
		struct SparseTraits
		{
			//!Vectors are usually dense.
			using Vector = Eigen::Matrix<double,Eigen::Dynamic,1>;
			//! Our problems are from /R^n to /R^n
			using IterationFunction= std::function <Vector (Vector const &)>;			
			using Matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
			
			//! Computes the distance between two vectors.
			static double distance (Vector const & current, Vector const & previous)
			{
				return (current-previous).norm() ;
				}
			
			//! Prints the vector \p v to the output stream \p OS
			static std::ostream & print (const Vector& v, std::ostream & OS)
			{
				for (std::size_t i = 0; i < v.size() ; ++i) OS << *(v.data()+i) << " ";
				OS << std::endl ;
				return OS;
			}
			
			//This one in case I want to use a sparse vector
			
			// static std::ostream & print (const Vector& v, std::ostream & OS)
			// {
			// for (Vector::InnerIterator it(v); it; ++it)
			// {
			// OS << it.index() <<" ";
			// OS << it.value() <<std::endl;
			// return OS;
			// }
			// }
			
		};
		
	}
	
	
	
	
#endif /* _FIXEDPOINTTRAITS_HPP_ */
