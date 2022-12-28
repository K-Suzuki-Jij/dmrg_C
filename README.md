# dmrg_C

Archived programs for calculating one-dimensional models about condensed matter physics using density matrix renormalization group (DMRG) and exact diagonalization.

## How to Use
First, you go enter `dmrg_C/lib` and do make

```
cd lib
make
```

If you use Linux systems, please comment out Mac OS clang commands 
https://github.com/K-Suzuki-Jij/dmrg_C/blob/afa53ab02f8b1cc591a355718f39988618e1f9f9/lib/makefile#L7-L9

and then uncomment Linux command.
https://github.com/K-Suzuki-Jij/dmrg_C/blob/afa53ab02f8b1cc591a355718f39988618e1f9f9/lib/makefile#L11-L13

We check the programes work on Ubuntu 20.04.5 LTS (GNU/Linux 5.15.0-56-generic x86_64)

Next, you go enter any directries under `dmrg_C/main/dmrg`.
Here, we go enter `dmrg_C/main/dmrg/xxz_vf` and do make

```
cd dmrg_C/main/dmrg/xxz_vf
make
```

Again here, If you use Linux systems, please comment out Mac OS clang commands 
https://github.com/K-Suzuki-Jij/dmrg_C/blob/afa53ab02f8b1cc591a355718f39988618e1f9f9/main/dmrg/xxz_vf/makefile#L14-L17

and then uncomment Linux command.
https://github.com/K-Suzuki-Jij/dmrg_C/blob/afa53ab02f8b1cc591a355718f39988618e1f9f9/main/dmrg/xxz_vf/makefile#L19-L22

Finally, you can execute program appeared in `dmrg_C/main/dmrg/xxz_vf/build/1D_XXZ_VF_DMRG.out`

```
./build/1D_XXZ_VF_DMRG.out
```

Calculation log will be appeared.
<img width="1385" alt="Screenshot 2022-12-29 at 2 00 15" src="https://user-images.githubusercontent.com/78338408/209846779-655dc791-7965-4853-881f-cd38356adac2.png">

Results will be written in `dmrg_C/main/dmrg/xxz_vf/result` and status will be written in `dmrg_C/main/dmrg/xxz_vf/SML_out`

## Programs by DMRG
* `dmrg_C/main/dmrg/hubbard_vf`: Hubbard model under the longitudinal fields.
* `dmrg_C/main/dmrg/klm_tvf`: Kondo lattice model under the transeverse fields.
* `dmrg_C/main/dmrg/klm_vf`: Kondo lattice model under the longitudinal fields.
* `dmrg_C/main/dmrg/tklm_vf`: Two-channel Kondo lattice model under the longitudinal fields.
* `dmrg_C/main/dmrg/xxz_vf`: XXZ model under the longitudinal fields.

## Programs by Exact Diagonalization
This is under construction. (maybe forever)
* `dmrg_C/main/exact_diag/klm_tvf`: Kondo lattice model under the transeverse fields.
* `dmrg_C/main/exact_diag/klm_vf`: Kondo lattice model under the longitudinal fields.

## Operation Check
* Mac OS 13.1 on Apple M1 Max
* Ubuntu 20.04.5 LTS on Ryzen 7950x


