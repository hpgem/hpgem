/* 
 * File:   MpiContainer.hpp
 * Author: brinkf
 *
 * Created on October 30, 2014, 4:31 PM
 */

#ifndef MPICONTAINER_HPP
#define	MPICONTAINER_HPP

class MpiContainer {
public:
    MpiContainer& Instance(){
        static MpiContainer theInstance();
        return theInstance;
    }
    
    int getProcessorID();
    int getCommSize();
#if HPGEM_USE_MPI
    MPI::Intracomm& getComm();
    
    template<class T>
    void broadcast(T&);
    
    template<class T>
    void send(T&);
#endif
    
private:
    MpiContainer();
    MpiContainer(const MpiContainer& orig)=delete;
    virtual ~MpiContainer();

};

#endif	/* MPICONTAINER_HPP */

