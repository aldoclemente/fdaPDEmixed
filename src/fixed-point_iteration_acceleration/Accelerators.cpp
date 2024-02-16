//! @file Accelerators.cpp
//! @brief Implementation of the acceleration algorithms
#include "Accelerators.hpp"
#include <Eigen/QR>

namespace FixedPoint
{	
	Traits::Vector  ASecantAccelerator::operator()(const std::deque < Vector > & past)
	{
		assert (!past.empty());
		Vector xNew (dimension);
		if (firstTime)
		{
			xNew = phi (past.back());
			phiOld = xNew;
			deltaXOld = xNew - past.back();
			firstTime = false;
		}
		else
		{
			Vector tmp (dimension), deltaXNew (dimension);
			xNew = phi (past.back());
			// compute \Delta x_n = phi(x_n)-x_n
			deltaXNew = xNew - past.back();
			
			// \Delta x_n - \Delta x_{n-1}
			tmp = deltaXNew - deltaXOld;
			// Get ||\Delta x_n - \Delta x_{n-1}||^2
			double norm2=tmp.squaredNorm();
			//to do: check norm2
			double factor = tmp.dot(deltaXNew);
			factor /= norm2;
			// Save phi(X_n) for later use
			tmp=xNew;
			// scale (phi(x_n)-\phi(x_{n-1}) factor and subtract
			// 
			for (std::size_t i=0;i<dimension;++i)
			xNew.coeffRef(i) -= factor * (xNew.coeffRef(i) - phiOld.coeffRef(i));
			phiOld=tmp;
			deltaXOld=deltaXNew;
		}
		
		return xNew ;
		
	}
	
	Traits::Vector  AndersonAccelerator::operator()(const std::deque < Vector > & past){
		
		assert (!past.empty());
		Vector f = phi(past.back()) - past.back() ;
		
		
		if ( memory == 1 ) //this case is called simple mixing
		{
			return (past.back() +  mixingParameter * f);
		}
		
		if ( past.size() == 1 )
		{
			fOld = f;
			return (past.back() +  mixingParameter * f);
		}
		
		else
		{
			X.push_back ( past[past.size()-1] - past[past.size()-2]) ;
			F.push_back ( f - fOld);
			fOld = f;
			//Calculating the solution using QR decomposition for solving least square problem
			//....Theoretical note about permutation matrix will be in the report....
			return (past.back() + mixingParameter*f - ( X.getMatrix() + mixingParameter*F.getMatrix())*F.getMatrix().colPivHouseholderQr().solve(f)) ;
		}
	}
	
}