# Notes on testing in hpGEM ##


## Unit tests ##

Unit test should, traditionally, test an independent unit of code and mock any 
outside parts that are needed. This however assumes that it is practically 
possible to mock outside dependencies. Such a clean separation is hard to 
achieve in an architecture with layered components. For example, to 
independently unit test the construction of a global matrix, one would need to
both mock a mesh and the associated basis functions. Mocking all these 
dependencies would require a large amount of work both to setup the mocks, not 
speaking of the maintenance required if interfaces change. The large cost in 
mocking defeats in our opinion the benefits from separating the components.

## Self Tests ##

The self tests are a mix of 'integration tests' and 'system tests'. The 
integration tests are nothing special, they test a larger part of the code or 
the interaction with dependencies. The 'system test's need to test the use of 
hpGEM.

Testing the whole system requires solving some kind of finite element problem 
and checking the output. The question with such an approach is, when is the 
numerical output correct? In a paper one would theoretically derive the 
convergence rate say h^(p+1) and then see if the errors in the numerical result
are close enough to this rate. For a system test this is not practical, as the 
criterion 'close enough' is not well defined enough for a compute to check. Thus
we take an alternative approach. We run the code once and verify the correct 
convergence rate manually. With the correctness verified the computed solutions
(with numerical errors included) are stored in the code and any subsequent run 
should give the exact same result.

If these system test fails there are therefore two possibilities:
  1. There is a bug.
  2. Something in the code changed so that a valid but different solution will
     be generated.
To distinguish between these two cases, one would need to look at the output and
compute the convergence rate. If both the output and the convergence rates look 
good, then the expected output may be updated, while if they are not it should 
be considered as a bug.
