module modtypes

    implicit none

    integer, parameter :: SP_ION = 0
    integer, parameter :: SP_ELECTRON = 1
    
    type externalTraits
        double precision :: temperature
        double precision, dimension(3) :: field(3)
    end type externalTraits

    type sampleTraits
        integer :: w, l, h
        integer :: nElectrons;
        double precision, dimension(3) :: easyAxis
    end type sampleTraits

    type interactionTraits
        double precision :: ee, eeCutOff
        double precision :: ei, eiCutOff
        double precision :: ii, iiCutOff
        double precision :: pauli
    end type interactionTraits

end module modtypes
