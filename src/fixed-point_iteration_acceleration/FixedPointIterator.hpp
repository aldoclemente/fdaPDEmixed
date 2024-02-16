//! @file FixedPointIterator.hpp
//! @brief Header file containing the class for solving a fixed-point problem

#ifndef SRC_FIXEDPOINTITERATOR_HPP_
	#define SRC_FIXEDPOINTITERATOR_HPP_
	
	#include <memory>
	#include <limits>
	#include "Accelerators.hpp"
	
	namespace FixedPoint
	{
		//! @brief Options to be composed into FixedPointIterator
		struct FixedPointOptions
		{
			//! Tolerance for the distance between two consecutive iterates to stop the iterations
			double tolerance=1.e-8;
			//! Maximum iterations
			unsigned int maxIter=100;
			//! number of data to keep stored in memory
			unsigned int memory = maxIter;
			//! number of data to print with FixedPointIterator::printHistory function
			unsigned int printMemory = 30;
		};
		
		//! @brief The core of the interface for solving a fixed-point problem
		//!
		//! Composes with a std::unique_ptr the other core class Iterator;
		//! contains several options and the history of the iteration process.
		class FixedPointIterator: public Traits
		{
			public:
			
			//TOFIX
			// This code has lead me to an Internal Compiler Error. Making the class default constructable seems good.
			// FixedPointIterator(std::unique_ptr <Iterator> && iter, const FixedPointOptions& opt=FixedPointOptions{}):
			// iterator(std::move(iter)), options{opt}, iteration{0}, history{} {}
			
			//! @brief Default constructor
			//!
			//!To have a meaningful object FixedPointIterator::setIterator will be called after this
			FixedPointIterator():iterator(nullptr), options{FixedPointOptions{}}, iteration{0}, history{} {}
			
			//! Gets #options (const version)
			FixedPointOptions getOptions() const {
				return options;
			}
			//! Gets #options
			FixedPointOptions & getOptions () {
				return options;
			}
			
			//! Sets #iterator
			void setIterator( std::unique_ptr<Iterator> && iter){iterator=std::move(iter);}
			//! Gets #iterator (const version)
			const Iterator & getIterator() const {return *iterator;}
			//! Gets #iterator
			Iterator & getIterator(){return *iterator;}
			//! Gets #iterator (pointer version)
			std::unique_ptr<Iterator> & getIterator_p(){return iterator;}

			
			
			//! @brief Calculates the sequence of iterates until convergence up to FixedPointOptions::tolerance or maximum iteration is reached
			//!
			//! If it starts from the first iteration, a default initial vector of zeros is used.
			//! @return True if converged
			bool compute ();
			//! @brief Calculates the sequence of iterates until convergence up to FixedPointOptions::tolerance or maximum iteration is reached
			//!
			//! @param x0 suggested initial value: it is used only if there are no past iterations.
			//! @return True if converged
			bool compute ( Vector const & x0 ) {if ( history.empty() ) history.emplace_back ( x0 ); return compute(); }
			
			//! @brief Resets the state of the object
			//!
			//! Resets the object pointed by #iterator, #history and #iteration 
			void reset() { 
				iterator->reset();
				history.clear();
				history.shrink_to_fit();
				iteration = 0;				
			}
			
			//! Computes last residual and prints useful infos
			void printResult ( std::ostream & OS = std::cout) const;
			
			//! @brief Prints #history to output stream \p OS, taking into consideration
			//! member FixedPointOptions::printMemory of #options
			void printHistory ( std::ostream & OS = std::cout) const {
				unsigned int m = std::min (history.size(), static_cast<long unsigned int>(options.printMemory));
				std::cout<< "Last " << m << " values:\n";
				for (int i = 0 ; i < m ; ++i) Traits::print ( history [i] , OS) ;
			}
			
			//! @brief Computes last residuals, taking into account #options
			//!
			//! It is a bit costly from the computational point of view, so use it with care
			void printResidualHistory ( std::ostream & OS = std::cout) const {
				unsigned int m = std::min (history.size(), static_cast<long unsigned int>(options.printMemory));
				std::cout<< "Last " << m << " residual norms:\n";
			for (int i = history.size()-m ; i < history.size() ; ++i) OS << distance( iterator->getIterationFunction()( history[i] ), history [i]) << std::endl ;
			}
			
			//! Gets #history
			const std::deque < Vector > & getHistory () const { return history ;}
			
			//! Gets #iteration (const version)
			const unsigned int getIteration () const {return iteration ;}
			
			//! Gets #iteration 
			unsigned int & getIteration () {return iteration ;}
			
			private:
			
			//! A std::unique_ptr to an Iterator object that produces the iterates
			std::unique_ptr<Iterator> iterator;
			//! Options
			FixedPointOptions options;
			//! @brief History of iterations
			//!
			//! The maximum number it can contain is established by FixedPointOptions::memory member
			std::deque < Vector > history;
			//! Number of computed iterations
			unsigned int iteration;
		};
		
		
	}
#endif /* SRC_FIXEDPOINTITERATOR_HPP_ */
