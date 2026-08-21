# Electron Density (Ne) Calculation and Update in MOCASSIN

This document provides a detailed trace of how the electron density (`Ne`) is calculated and where its iteration value is updated within the codebase. It also outlines the potential effects of dust grains on the overall simulation.

## 1. Primary Calculation of Electron Density (`Ne`)

The fundamental calculation and updating of the electron density takes place in the `ionBalance` subroutine, which handles the iterative ionization balance for the gas in a given cell.

**File:** `source/update_mod.f90`
**Subroutine:** `ionBalance` (and previously recursive `ionBalance`)
**Variables Involved:**
- `NeUsed`: A temporary variable accumulating the electron density for the current cell before assigning it to the grid arrays.
- `grid%elemAbun`: Element abundance array.
- `ionDenUsed`: The fractional abundance of a particular ion.
- `grid%Hden(cellP)`: The total hydrogen density in the cell.
- `grid%Ne(cellP)`: The grid array holding the electron density.
- `NeTemp(cellP)`: A temporary electron density array used to assist with MPI implementation and next iterations.

### Calculation Steps

The calculation sums the free electrons contributed by all elements and their respective ionization stages. For a given element and ionization stage `ion` (where `ion=1` is neutral, `ion=2` is singly ionized, etc.), the number of free electrons contributed is `(ion-1)`. This is multiplied by the elemental abundance and the ionic fraction, and finally scaled by the total hydrogen density.

Here is the exact code snippet showing how `NeUsed` is calculated and updated:

```fortran
            ! calculate new Ne
          NeUsed = 0.
            do elem = 1, nElements
                do ion = 2, min(elem+1, nstages)
                    if (lgElementOn(elem)) then
                      if( ionDenUsed(elementXref(elem),ion) >= 1.e-10) &
                           & NeUsed = NeUsed + (ion-1)*&
                           &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                           &ionDenUsed(elementXref(elem), ion)
                    end if

                end do
            end do

            NeUsed = NeUsed * grid%Hden(cellP)
```

### Iteration Value Update

Immediately following this calculation, if `NeUsed` falls to zero, it is bumped up to a minimum value of `1.0`. Then, the newly computed value is saved back into the `grid` struct and a temporary `NeTemp` array which is used for the next iteration and MPI communication.

```fortran
            if (NeUsed == 0.) then
               NeUsed = 1.
            end if

            if (LgNeInput) then
               correction = NeUsed/grid%NeInput(cellP)
               grid%Hden(cellP) = grid%Hden(cellP)/correction
            end if

            ! this was added to help MPI implementation
            NeTemp(cellP) = NeUsed

            grid%Ne(cellP) = NeUsed
```

*(Note: There is also an iterative version of this block commented out in the same file `source/update_mod.f90`, which exhibits the same `NeUsed` updating logic.)*

---

## 2. Calculation of Heavy Metal Free Electron Density for Opacity

Another calculation related to electron density happens when the code needs to evaluate the free-free (bremsstrahlung) opacity. Here, the contribution specifically from heavy metals is summed up.

**File:** `source/ionization_mod.f90`
**Subroutine:** `eDenSum`
**Variables Involved:**
- `eDenFFSum`: A global/module variable accumulating the heavy element free electron density.
- `density(n, i+1)`: The number density of the specific ion.

### Calculation Steps

The subroutine iterates over all elements starting from `n=3` (skipping Hydrogen and Helium) and their respective ionization stages to compute `eDenFFSum`. The formula used includes a factor of `i*i` where `i` represents the charge of the ion, effectively weighting the density to compute the parameter used for the gaunt factor in free-free emission/absorption.

```fortran
    ! this subroutine sums up the free electron density over all species
    subroutine eDenSum()
        implicit none

        ! local variables
        integer :: n, i                                           ! counters

        ! sum the heavy metal free electron density
        eDenFFSum = 0.
        do n = 3, nElements
            do i = 1, min(n, nstages-1)
                eDenFFSum = eDenFFSum + float(i*i)*density(n, i+1)
            end do
        end do
    end subroutine eDenSum
```

This `eDenFFSum` is subsequently used inside the `addOpacity` subroutine to evaluate the free-free opacity factor:

```fortran
        ! hydrogen helium and heavy element brems (free-free) opacity,
        ! assuming hydrogen ff gaunt factors
        ! xSecArray is missing factor of 1.e-20 to avoid underflow.
        fac1 = (NeUsed/1.e10) * ( (density(1,2) + density(2, 2) + &
             & 4.*density(2, 3) + eDenFFSum)/1.e10 ) / sqrTeUsed
```

---

## 3. Potential Contributions and Effects from Dust Grains

In MOCASSIN, the presence of dust is controlled by the `lgDust` logical flag. While dust grains do not directly emit "free electrons" that modify the `NeUsed` variable calculated in `ionBalance`, they significantly alter the radiation field, thermal balance, and therefore implicitly affect the ionization balance (and the final electron density).

When `lgDust` is `.true.`, several key mechanisms take effect:

1. **Dust Opacity (Absorption and Scattering):**
   Dust heavily influences the transfer of radiation. The opacities are calculated in `source/dust_mod.f90` inside `dustOpacity()`.
   ```fortran
           grid%scaOpac(grid%active(iP,jP,kP),i) = grid%Ndust(grid%active(iP,jP,kP))*&
                & xSecArray(dustScaXsecP(0,1)+i-1)
           grid%absOpac(grid%active(iP,jP,kP),i) = grid%Ndust(grid%active(iP,jP,kP))*&
                & xSecArray(dustAbsXsecP(0,1)+i-1)
   ```
   By absorbing ionizing photons, dust decreases the available radiation for gas photoionization, indirectly lowering the fractional ion abundances (`ionDenUsed`) and thereby lowering the resulting electron density `NeUsed`.

2. **Photoelectric Heating:**
   Dust grains can eject electrons when struck by high-energy photons. While the number of these electrons is negligible compared to gas photoionization (hence not added to `NeUsed`), the kinetic energy they carry is a major heating source for the gas. This is computed in `thermBalance()` inside `source/update_mod.f90`:
   ```fortran
            if (lgDust .and. lgPhotoelectric) then
               coolInt = coolInt+gasDustColl_g
               heatInt = heatInt+photoelHeat_g
            end if
   ```
   This alters the gas temperature (`TeUsed`). Since recombination rates (`alphaTot`) strongly depend on the temperature, altering the temperature shifts the ionization balance, once again indirectly affecting `NeUsed`.

3. **Dust-Gas Collisions (Heating/Cooling):**
   Collisions between free electrons/ions and dust grains exchange energy, which also acts as a heating/cooling term for the gas (`gasDustColl_g`). This provides an additional pathway through which dust modifies the electron temperature, changing the rates that dictate the electron density.

In summary, dust fundamentally impacts the environmental conditions (radiation field and gas temperature) that dictate the iterative calculation of the electron density, even though it does not explicitly add a direct density term to the `NeUsed` variable.
