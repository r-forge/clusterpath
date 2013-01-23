library(inline)
funx <- cfunction(signature(), '
                /* creating a pointer to a vector<int> */
                std::vector<int>* v = new std::vector<int> ;
                v->push_back( 1 ) ;
                v->push_back( 2 ) ;
                
                /* wrap the pointer as an external pointer */
                /* this automatically protected the external pointer from R garbage 
                   collection until p goes out of scope. */
                Rcpp::XPtr< std::vector<int> > p(v, true) ;
                
                /* return it back to R, since p goes out of scope after the return 
                   the external pointer is no more protected by p, but it gets 
                   protected by being on the R side */
                return( p ) ;
        ', Rcpp=TRUE, verbose=FALSE)
        xp <- funx()
