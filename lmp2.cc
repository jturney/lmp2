/*
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libdpd/dpd.h>
#include "psifiles.h"
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/mints.h>
#include <libmints/view.h>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

INIT_PLUGIN

namespace psi{ namespace lmp2{

void denom(const int& nirrep,
           const Dimension& occOrbsPI,
           const Dimension& virOrbsPI,
           const Dimension& occOffset,
           const Dimension& virOffset);

void pair_energies(IntegralTransform& ints, const int& nirrep, const Dimension& occOrbsPI, double** epair);
void print_pair_energies(const int& nirrep, const Dimension& occOrbsPI, double* emp2);

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "LMP2" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" PsiReturnType
lmp2(Options &options)
{
    fprintf(outfile,"\n\n *******************************************************************************\n");
    fprintf(outfile,    " *                                  Local MP2                                  *\n");
    fprintf(outfile,    " *                    by Justin Turney and Brandon Magers                      *\n");
    fprintf(outfile,    " *           Parts taken from ccsort and ccenergy by Daniel Crawford           *\n");
    fprintf(outfile,    " *******************************************************************************\n\n");

    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");

    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");

    // Quickly check that there are no open shell orbitals here...
    int h, i, j, a, b;
    int numAOcc, numBOcc, numAVir, numBVir;
    double eSCF;
    double *occEvals, *virEvals;

    int nirrep          = Process::environment.wavefunction()->nirrep();

    // The localization process removes symmetry, make sure the user knows that.
    if (nirrep != 1) {
        fprintf(outfile, "    Localization removes symmetry, this run is not being conducted in C1!\n");
        return Failure;
    }
    SharedVector aEvals = Process::environment.wavefunction()->epsilon_a();
    SharedVector bEvals = Process::environment.wavefunction()->epsilon_b();
    char **labels       = Process::environment.molecule()->irrep_labels();
    int nmo             = Process::environment.wavefunction()->nmo();
    Dimension mopi      = Process::environment.wavefunction()->nmopi();
    Dimension clsdpi    = Process::environment.wavefunction()->doccpi();
    Dimension openpi    = Process::environment.wavefunction()->soccpi();
    Dimension frzcpi    = Process::environment.wavefunction()->frzcpi();
    Dimension frzvpi    = Process::environment.wavefunction()->frzvpi();
    Dimension occOrbsPI(nirrep);
    Dimension virOrbsPI(nirrep);
    Dimension occOffset(nirrep);
    Dimension virOffset(nirrep);
    numAOcc = 0; numBOcc = 0; numAVir = 0; numBVir = 0;
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    for(h = 0; h < nirrep; ++h){
        occOrbsPI[h] = clsdpi[h] - frzcpi[h];
        virOrbsPI[h] = mopi[h]   - clsdpi[h] - frzvpi[h];
        numAOcc += occOrbsPI[h];
        numAVir += virOrbsPI[h];
    }

    // form offsets
    int ocount = occOrbsPI[0], vcount = virOrbsPI[0];
    for (h=1; h<nirrep; h++) {
        occOffset[h] = ocount;
        ocount += occOrbsPI[h];
        virOffset[h] = vcount;
        vcount += virOrbsPI[h];
    }

    occEvals = new double[numAOcc];
    virEvals = new double[numAVir];
    fprintf(outfile, "\n\n\tIrrep  Core  Docc  Socc  Occ  Vir\n");
    fprintf(outfile,     "\t===============================================\n");
    for(h = 0; h < nirrep; ++h){
       fprintf(outfile, "\t %3s   %3d   %3d   %3d   %3d  %3d\n",
                             labels[h], frzcpi[h], clsdpi[h], openpi[h],
                             occOrbsPI[h], virOrbsPI[h]);
    }
    fprintf(outfile,     "\t===============================================\n\n");
    aOccCount = 0; aVirCount = 0;
    for(h = 0; h < nirrep; ++h){
        for(a = frzcpi[h]; a < clsdpi[h]; ++a) occEvals[aOccCount++] = aEvals->get(h, a);
        for(a = clsdpi[h]; a < mopi[h]; ++a) virEvals[aVirCount++] = aEvals->get(h, a);
    }
    if(print > 2){
        for(i = 0; i < numAOcc; ++i)
            fprintf(outfile, "\toccEvals[%2d] = %10.6f\n", i, occEvals[i]);
        fprintf(outfile, "\n");
        for(i = 0; i < numAVir; ++i)
            fprintf(outfile, "\tvirEvals[%2d] = %10.6f\n", i, virEvals[i]);
        fprintf(outfile, "\n");
    }

    // => Two-electron integral transformation <=

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    IntegralTransform ints(wfn, spaces);

    {
        // Use the IntegralTransform object's DPD instance, for convenience
        dpd_set_default(ints.get_dpd_id());

        // Since we're doing multiple separate transformations, keep the SO ints around:
        ints.set_keep_dpd_so_ints(1);

#if 0
        // Generate the integrals in various spaces in chemists' notation
        ints.set_dpd_int_file(PSIF_CC_AINTS);
        ints.set_aa_int_name("A (ik|jl)");
        fprintf(outfile, "\tTransforming integrals into types of A (ik|jl) ...\n\n"); fflush(outfile);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
        fprintf(outfile, "\n");

        ints.set_dpd_int_file(PSIF_CC_CINTS);
        ints.set_aa_int_name("C (ij|ab)");
        fprintf(outfile, "\tTransforming integrals into types of C (ij|ab) ...\n\n"); fflush(outfile);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
        fprintf(outfile, "\n");

        // This is the last transformation, we can delete the SO ints
        ints.set_keep_dpd_so_ints(0);
        ints.set_dpd_int_file(PSIF_CC_DINTS);
        ints.set_aa_int_name("D (ia|jb)");
        fprintf(outfile, "\tTransforming integrals into types of D (ia|jb) ...\n\n"); fflush(outfile);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

        fprintf(outfile, "\n\tIntegral transformations complete.\n\n");
#endif

        // Open files we need
        psio->open(PSIF_CC_OEI, PSIO_OPEN_NEW);
        psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
        psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
        psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
        psio->open(PSIF_CC_DENOM, PSIO_OPEN_OLD);
        psio->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);

        // => Re-sort the chemists' notation integrals to physicists' notation (pq|rs) = <pr|qs> <=
        dpdbuf4 A, C, D;
#if 0
        // Form A <ij|kl> from A (ik|jl)
        dpd_buf4_init(&A, PSIF_CC_AINTS, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>=O]+"), ID("[O>=O]+"), 0, "A (ik|jl)");
        fprintf(outfile, "\tSorting A (ik|jl) to A <ij|kl> ... "); fflush(outfile);
        dpd_buf4_sort(&A, PSIF_CC_AINTS, prqs, ID("[O,O]"), ID("[O,O]"), "A <ij|kl>");
        fprintf(outfile, "done\n"); fflush(outfile);
        dpd_buf4_close(&A);

        // Form C <ia|jb> from C (ij|ab)
        dpd_buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "C (ij|ab)");
        fprintf(outfile, "\tSorting C (ij|ab) to C <ia|jb> ... "); fflush(outfile);
        dpd_buf4_sort(&C, PSIF_CC_CINTS, prqs, ID("[O,V]"), ID("[O,V]"), "C <ia|jb");
        fprintf(outfile, "done\n"); fflush(outfile);
        dpd_buf4_close(&C);

        // Form D <ij|ab> from D (ia|jb)
        dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "D (ia|jb)");
        fprintf(outfile, "\tSorting D (ia|jb) to D <ij|ab> ... "); fflush(outfile);
        dpd_buf4_sort(&D, PSIF_CC_DINTS, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
        fprintf(outfile, "done\n");
        dpd_buf4_close(&D);
#endif

        // Form D 2<ij|ab> - <ij|ba>
        dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
        fprintf(outfile, "\tForming D 2<ij|ab> - <ij|ba> ... "); fflush(outfile);
        dpd_buf4_copy(&D, PSIF_CC_DINTS, "D 2<ij|ab> - <ij|ba>");
        fprintf(outfile, "<ij|ab> done ... "); fflush(outfile);
        {
            dpdbuf4 Dba;
            dpd_buf4_init(&Dba, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1, "D <ij|ab>");
            dpd_buf4_axpy(&Dba, &D, 1.0);
            fprintf(outfile, "<ij|ab> - <ij|ba> done ... "); fflush(outfile);
            dpd_buf4_close(&Dba);
        }
        fprintf(outfile, "done\n");
        dpd_buf4_close(&D);
    }

    // Form the Fock matrices
    {
        fprintf(outfile, "\n\tViewing Fock matrices ... "); fflush(outfile);

        // Looks like libtrans did this for us. Load it up.
        SharedMatrix f(new Matrix("f(m,n)", mopi, mopi));
        f->load(psio, PSIF_OEI, PSIF_MO_FOCK, nmo);

        // View out our frozen orbitals
        View view_oei_ij(f, occOrbsPI, occOrbsPI, frzcpi, frzcpi);
        SharedMatrix fij = view_oei_ij();
        fij->set_name("fIJ");

        View view_oei_ab(f, virOrbsPI, virOrbsPI, frzcpi + occOrbsPI, frzcpi + occOrbsPI);
        SharedMatrix fab = view_oei_ab();
        fab->set_name("fAB");

        dpdfile2 F;
        dpd_file2_init(&F, PSIF_CC_OEI, 0, ID('O'), ID('O'), "fIJ");
        fij->write_to_dpdfile2(&F);
        dpd_file2_close(&F);

        dpd_file2_init(&F, PSIF_CC_OEI, 0, ID('V'), ID('V'), "fAB");
        fab->write_to_dpdfile2(&F);
        dpd_file2_close(&F);

        dpd_file2_init(&F, PSIF_CC_OEI, 0, ID('O'), ID('O'), "fij");
        fij->write_to_dpdfile2(&F);
        dpd_file2_close(&F);

        dpd_file2_init(&F, PSIF_CC_OEI, 0, ID('V'), ID('V'), "fab");
        fab->write_to_dpdfile2(&F);
        dpd_file2_close(&F);

        fprintf(outfile, "done\n\n");
    }

    // => Form denominators <=
    {
        fprintf(outfile, "\tForming denominators ... "); fflush(outfile);
        denom(nirrep, occOrbsPI, virOrbsPI, occOffset, virOffset);
        fprintf(outfile, "done\n\n"); fflush(outfile);
    }

    // => Form initial LMP2 amplitudes <=
    {
        dpdbuf4 D;
        dpdbuf4 T2;

        dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
        dpd_buf4_copy(&D, PSIF_CC_TAMPS, "LMP2 tIjAb");
        dpd_buf4_close(&D);

        dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
        dpd_buf4_init(&D, PSIF_CC_DENOM, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "dIjAb");
        dpd_buf4_dirprd(&D, &T2);
        dpd_buf4_close(&D);
        dpd_buf4_close(&T2);
    }

    // => Compute the LMP2 energy <=
    double energy = 0.0;
    double rms = 0.0;
    int conv = 0;
    int lmp2_maxiter = 1000;
    {
        dpdbuf4 D, T2, newT2;
        dpdfile2 fij, fab;

        dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D 2<ij|ab> - <ij|ba>");
        dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
        energy = dpd_buf4_dot(&D, &T2);
        dpd_buf4_close(&T2);
        dpd_buf4_close(&D);

        fprintf(outfile, "\tComputing LMP2 amplitudes:\n");
        fprintf(outfile, "\t==========================\n");
        fprintf(outfile, "\titer = %3d LMP2 Energy = %20.14lf\n", 0, energy); fflush(outfile);

        for (int iter=1; iter< lmp2_maxiter; iter++) {
            dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
            dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New LMP2 tIjAb Increment");
            dpd_buf4_close(&D);

            dpd_buf4_init(&newT2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb Increment");
            dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
            dpd_file2_init(&fij, PSIF_CC_OEI, 0, ID('O'), ID('O'), "fIJ");
            dpd_contract424(&T2, &fij, &newT2, 1, 0, 1, -1.0, 1.0);
            dpd_contract244(&fij, &T2, &newT2, 0, 0, 0, -1.0, 1.0);
            dpd_file2_close(&fij);

            dpd_file2_init(&fab, PSIF_CC_OEI, 0, ID('V'), ID('V'), "fAB");
            dpd_contract244(&fab, &T2, &newT2, 1, 2, 1, 1.0, 1.0);
            dpd_contract424(&T2, &fab, &newT2, 3, 1, 0, 1.0, 1.0);
            dpd_buf4_copy(&T2, PSIF_CC_TAMPS, "New LMP2 tIjAb");
            dpd_buf4_close(&T2);

            dpd_buf4_init(&D, PSIF_CC_DENOM, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "dIjAb");
            dpd_buf4_dirprd(&D, &newT2);
            dpd_buf4_close(&D);
            dpd_buf4_close(&newT2);

            dpd_buf4_init(&newT2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb");
            dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb Increment");
            dpd_buf4_axpy(&T2, &newT2, 1.0);
            dpd_buf4_close(&T2);

            dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D 2<ij|ab> - <ij|ba>");
            energy = dpd_buf4_dot(&D, &newT2);
            dpd_buf4_close(&D);

            dpd_buf4_close(&newT2);

            /* Check for convergence */
            dpd_buf4_init(&newT2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb");
            dpd_buf4_mat_irrep_init(&newT2, 0);
            dpd_buf4_mat_irrep_rd(&newT2, 0);

            dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
            dpd_buf4_mat_irrep_init(&T2, 0);
            dpd_buf4_mat_irrep_rd(&T2, 0);

            rms = 0.0;
            for(int row=0; row < T2.params->rowtot[0]; row++)
              for(int col=0; col < T2.params->coltot[0]; col++)
                rms += (newT2.matrix[0][row][col] - T2.matrix[0][row][col]) *
                  (newT2.matrix[0][row][col] - T2.matrix[0][row][col]);

            dpd_buf4_mat_irrep_close(&T2, 0);
            dpd_buf4_mat_irrep_close(&newT2, 0);
            dpd_buf4_close(&T2);
            dpd_buf4_close(&newT2);

            rms = sqrt(rms);

            fprintf(outfile, "\titer = %3d LMP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms); fflush(outfile);

            double convergence = 1.0e-8;
            if(rms < convergence) {
              conv = 1;
              fprintf(outfile, "\n\tLMP2 Iterations converged.\n");
              break;
            }
            else {
              dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb");
              dpd_buf4_copy(&T2, PSIF_CC_TAMPS, "LMP2 tIjAb");
              dpd_buf4_close(&T2);
            }
        }

        if(!conv) {
          fprintf(outfile, "\n\tLMP2 Iterative procedure failed.\n");
          throw ConvergenceError<int>("LMP2 interative procedure failed.", lmp2_maxiter, 1.0e-8, rms, __FILE__, __LINE__);
        }

        fprintf(outfile, "\tLMP2 Correlation Energy = %20.14f\n", energy);
        fprintf(outfile, "\tLMP2 Total Energy       = %20.14f\n\n", energy+Process::environment.globals["HF TOTAL ENERGY"]);
        fflush(outfile);

        dpd_buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
        dpd_buf4_close(&T2);

        double* epair;
        pair_energies(ints, nirrep, occOrbsPI, &epair);
        print_pair_energies(nirrep, occOrbsPI, epair);
    }

    // Close files we opened.
    psio->close(PSIF_CC_OEI, 0);
    psio->close(PSIF_CC_AINTS, 0);
    psio->close(PSIF_CC_CINTS, 0);
    psio->close(PSIF_CC_DINTS, 0);
    psio->close(PSIF_CC_DENOM, 0);
    psio->close(PSIF_CC_TAMPS, 0);

    return Success;
}

}} // End Namespaces
