def write_plumed_file(p0, p1, protein_IDs, lig_IDs, extent=0.60, extent_buffer=0.15,
                      l_proj=0.5, u_proj=4.0, beta_cent=1.5, s_cent=2,
                      deposition_pace=1000, print_pace=1000, write_ProjectionOnAxis=False):
    """
    Writes a standard wt fun-metaD plumed.dat file with units in Ångström and kcal/mol.
    """

    version = 1.0
    nm_to_A = 10.0
    kj_to_kcal = 0.239005736

    p0_str = ','.join(str(i) for i in p0)
    p1_str = ','.join(str(i) for i in p1)
    protein_str = f"{protein_IDs[0]}-{protein_IDs[-1]}"
    lig_str = f"{lig_IDs[0]}-{lig_IDs[-1]}"

    with open('plumed.dat', 'w') as FILE:
        FILE.write('####################################\n')
        FILE.write('# plumed.dat for Funnel Metadynamics\n')
        FILE.write(f'# Written on {datetime.datetime.now()}\n')
        FILE.write(f'# By funnel_maker {version}\n')
        FILE.write('####################################\n')
        FILE.write('UNITS LENGTH=A ENERGY=kcal/mol\n')
        FILE.write('RESTART\n\n')

        if write_ProjectionOnAxis:
            FILE.write('LOAD FILE=ProjectionOnAxis.cpp\n\n')

        FILE.write(f'WHOLEMOLECULES STRIDE=1 ENTITY0={protein_str} ENTITY1={lig_str}\n\n')
        FILE.write(f'lig: COM ATOMS={lig_str}\n')
        FILE.write(f'p0: COM ATOMS={p0_str}\n')
        FILE.write(f'p1: COM ATOMS={p1_str}\n\n')

        FILE.write('pp: PROJECTION_ON_AXIS AXIS_ATOMS=p0,p1 ATOM=lig\n\n')

        FILE.write(f's_cent: CONSTANT VALUES={s_cent * nm_to_A:.1f}\n')
        FILE.write(f'beta_cent: CONSTANT VALUES={beta_cent:.1f}\n')
        FILE.write(f'wall_width: CONSTANT VALUES={extent * nm_to_A:.2f}\n')
        FILE.write(f'wall_buffer: CONSTANT VALUES={extent_buffer * nm_to_A:.2f}\n')
        FILE.write(f'lwall: LOWER_WALLS ARG=pp.proj AT={l_proj * nm_to_A:.1f} '
                   f'KAPPA={20000.0 * kj_to_kcal:.1f} EXP=2 EPS=1\n')
        FILE.write(f'uwall: UPPER_WALLS ARG=pp.proj AT={u_proj * nm_to_A:.1f} '
                   f'KAPPA={20000.0 * kj_to_kcal:.1f} EXP=2 EPS=1\n\n')

        FILE.write('MATHEVAL ...\n')
        FILE.write('  LABEL=wall_center\n')
        FILE.write('  ARG=pp.proj,s_cent,beta_cent,wall_width,wall_buffer\n')
        FILE.write('  VAR=s,sc,b,h,f\n')
        FILE.write('  FUNC=h*(1./(1.+exp(b*(s-sc))))+f\n')
        FILE.write('  PERIODIC=NO\n')
        FILE.write('... MATHEVAL\n\n')

        FILE.write('scaling: CONSTANT VALUES=1.0\n')
        FILE.write(f'spring: CONSTANT VALUES={1000.0 * kj_to_kcal:.1f}\n\n')

        FILE.write('MATHEVAL ...\n')
        FILE.write('  LABEL=wall_bias\n')
        FILE.write('  ARG=pp.ext,spring,wall_center,scaling\n')
        FILE.write('  VAR=z,k,zc,sf\n')
        FILE.write('  FUNC=step(z-zc)*k*(z-zc)*(z-zc)/(sf*sf)\n')
        FILE.write('  PERIODIC=NO\n')
        FILE.write('... MATHEVAL\n\n')

        FILE.write('finalbias: BIASVALUE ARG=wall_bias\n\n')

        FILE.write('METAD ...\n')
        FILE.write('  LABEL=meta ARG=pp.proj,pp.ext\n')
        FILE.write(f'  SIGMA={0.025 * nm_to_A:.3f},{0.03 * nm_to_A:.3f} '
                   f'HEIGHT={1.5 * kj_to_kcal:.3f}\n')
        FILE.write(f'  PACE={deposition_pace} FILE=HILLS\n')
        FILE.write(f'  GRID_MIN={(l_proj - 0.5) * nm_to_A:.1f},0.0 '
                   f'GRID_MAX={(u_proj + 0.5) * nm_to_A:.1f},'
                   f'{(extent + extent_buffer + 0.2) * nm_to_A:.1f} '
                   'GRID_SPACING=0.05,0.6\n')
        FILE.write('  BIASFACTOR=10.0 TEMP=298\n')
        FILE.write('... METAD\n\n')

        FILE.write(f'PRINT ARG=* STRIDE={print_pace} FILE=COLVAR FMT=%8.4f\n')

