# 38: Silicon — Optimized Projection Functions Method (OPFM)

- Outline: *Construct optimal projection functions for Si by
    codiagonalizing the Bloch-space overlap matrices, then
    inspect the resulting codiagonalization matrix.*

- Directory: [`tutorial/tutorial38/`](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial38)

- Input Files

    - `silicon.scf` *The `pwscf` input file for ground
        state calculation*

    - `silicon.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `silicon.pw2wan` *Input file for `pw2wannier90`*

    - `silicon.win` *The `wannier90` input file*

!!! note
    OPFM requires `wannier90` to be built and linked against
    [`libcodiag`](https://github.com/wannier-developers/libcodiag).
    See `README.install` for the `CODIAG_ROOT`/`WANNIER90_CODIAG` (CMake)
    or `CODIAG_INCLUDE`/`CODIAG_LIB` (`make.inc`) build options. Without
    it, `wannier90` still builds and runs normally; setting `opfm = true`
    in the input file fails with an input error at run time.

1. Run `pwscf` to obtain the ground state of silicon

    ```bash title="Terminal"
    pw.x < silicon.scf > scf.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid.

    ```bash title="Terminal"
    pw.x < silicon.nscf > nscf.out
    ```

3. Run `wannier90` to generate a list of the required overlaps (written
    into the `silicon.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp silicon
    ```

4. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silicon.mmn`
    and `silicon.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < silicon.pw2wan > pw2wan.out
    ```

    `silicon.win` defines several candidate sets of trial projections in
    its `projections` block, and selects one at a time via
    `select_projections`. Before enabling OPFM, run plain
    Marzari-Vanderbilt (MV) minimization on three of them to see how
    sensitive the result is to the choice of trial orbitals.

    **(a)** The default selection, `select_projections 1 2 3 4`, uses
    four bond-centred s-orbitals. Run `wannier90` as-is

    ```bash title="Terminal"
    wannier90.x silicon
    ```

    Inspect `silicon.wout`: MV minimization converges cleanly to a
    low final spread.

    **(b)** Comment out the bond-centred selection and uncomment the
    correctly-oriented atom-centred sp³ orbitals instead

    ```vi title="Input file"
    !! select_projections 1 2 3 4         !! bond-centered s-orbitals
    !! select_projections 5-8             !! sp3 incorrect orientation
    select_projections 9-12            !! sp3 correct orientation
    ```

    Re-run `wannier90.x silicon` and inspect `silicon.wout` again.
    This converges to the same low spread as (a), reached from a
    different starting guess.

    **(c)** Now select the sp³ orbitals with the *incorrect* orientation
    instead

    ```vi title="Input file"
    !! select_projections 1 2 3 4         !! bond-centered s-orbitals
    select_projections 5-8             !! sp3 incorrect orientation
    !! select_projections 9-12            !! sp3 correct orientation
    ```

    Re-run `wannier90.x silicon` and compare `silicon.wout` to the
    previous two runs. This time MV minimization gets stuck in a
    local minimum, converging to a noticeably higher final spread
    despite starting from projections of the same symmetry type as
    (b), just misoriented relative to the bond network.

5. Leaving `select_projections` set to the incorrect sp³ orientation
    (5-8) from step 4c, enable OPFM by uncommenting the following line
    in `silicon.win`

    ```vi title="Input file"
    opfm = .true.
    ```

6. Run `wannier90` again

    ```bash title="Terminal"
    wannier90.x silicon
    ```

    With OPFM constructing the initial guess instead of a direct
    projection, the same misoriented sp³ selection now converges to the
    same low spread as the correctly-oriented case from step 4b. OPFM
    recovers from the bad orientation that trapped plain MV minimization
    in step 4c.

7. Switch to a redundant, over-complete set of trial projections: the
    s and p orbitals on Si1 and its four neighbors (20 trial projections
    for 4 Wannier functions). Also uncomment `opfm_write_w_matrix` to
    write out the codiagonalization matrix for plotting in the next
    step — it defaults to `false`, since it's an extra output most
    OPFM runs don't need.

    ```vi title="Input file"
    !! select_projections 5-8             !! sp3 incorrect orientation
    select_projections 13-16 , 21-36   !! s,p orbitals on Si1 and its 4 neighbors
    opfm_write_w_matrix = .true.
    ```

    Run `wannier90.x silicon` again. Here OPFM has real freedom to
    combine the 20 candidates into the best 4: compare the spread in the
    "Initial State" block of `silicon.wout` (right after OPFM
    construction, before any MV iterations) to the "Final State" block —
    OPFM alone already lands close to the same minimum found in every
    other successful case above, with full MV minimization only
    improving it slightly further.

    This run writes `silicon_opfm_w.mat`, the $M \times N$
    codiagonalization matrix $W$ mapping the $M=20$ raw trial
    projections onto the $N=4$ optimized projection functions.

8. Plot the codiagonalization matrix using the `plot.py` script
    provided in the tutorial directory. Its dependencies (`numpy` and
    `matplotlib`) are declared inline via
    [PEP 723](https://peps.python.org/pep-0723/) script metadata, so
    running it with [`uv run`](https://docs.astral.sh/uv/) — or simply
    executing it directly, since its shebang already invokes `uv run` —
    fetches them into an ephemeral environment automatically, with no
    separate `pip install` step needed. Any Python environment with
    `numpy`/`matplotlib` already installed works too.

    ```bash title="Terminal"
    ./plot.py w-matrix silicon_opfm_w.mat
    ```

    This produces `silicon_opfm_w.png`, a heatmap of $|[W^T]_{ij}|$
    (Projection Index on the x-axis, OPF Index on the y-axis) showing
    how each of
    the 20 raw trial projections is redistributed across the 4 optimized
    projection functions. Each row sums to 1 in $|W_{ij}|^2$ by
    construction, so the colour scale is fixed to $[0,1]$.

## Further ideas

- Compare the Wannier spreads with and without `opfm = true` to see how
    much the optimized projections improve on the raw trial projections.

## References

- [10.1103/PhysRevB.92.165134](https://doi.org/10.1103/PhysRevB.92.165134)
- [10.1103/PhysRevB.94.125151](https://doi.org/10.1103/PhysRevB.94.125151)
