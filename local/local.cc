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

#include <boost/algorithm/string.hpp>
#include <utility>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>

#include <libmints/mints.h>
#include "local.h"
#include <physconst.h>

using namespace boost;
using namespace psi;

namespace psi { namespace local {

    Local::Local(boost::shared_ptr<Wavefunction> wavefunction) :
        wavefunction_(wavefunction), basisset_(wavefunction->basisset()), C_USO_(wavefunction->Ca()), print_(0), debug_(0)
    {
        common_init();
    }
    Local::Local(boost::shared_ptr<Wavefunction> wavefunction, boost::shared_ptr<BasisSet> auxiliary) :
        wavefunction_(wavefunction), basisset_(wavefunction->basisset()), auxiliary_(auxiliary), C_USO_(wavefunction->Ca()), print_(0), debug_(0)
    {
        common_init();
    }
    Local::~Local()
    {
    }
    void Local::common_init()
    {
        boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));

        if (C_USO_->nirrep() > 1) {
            boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
            AO2USO_ = pet->aotoso();
        }

        C_AO_ = C_AO();
        L_AO_ = C_AO_;

        int nso = basisset_->nbf();
        S_ = SharedMatrix(new Matrix("S",nso,nso));
        X_ = SharedMatrix(new Matrix("S^+1/2",nso,nso));
        boost::shared_ptr<OneBodyAOInt> Sint(integral->ao_overlap());
        Sint->compute(S_);
        X_->copy(S_);
        X_->power(+1.0/2.0);
    }
    void Local::print(FILE* out)
    {
        fprintf(out, "  ==> Localization <==\n\n");

        basisset_->print_by_level(out,3);
        if (auxiliary_.get())
            auxiliary_->print_by_level(out,3);

        C_USO_->print(out);
        C_AO_->print(out);
        L_AO_->print(out);

        fprintf(out, "  => Metrics <=\n\n");
        fprintf(out, "  Boys:                %11.3E\n", boys_metric());
        fprintf(out, "  Edmiston-Ruedenberg: %11.3E\n", er_metric());
        fprintf(out, "  Pipek-Mezey:         %11.3E\n\n", pm_metric());

        if (Q_.get())
            Q_->print(out);

        if (domains_.size()) {
            fprintf(out, "  => Primary Domains <=\n\n");
            for (int i = 0; i < L_AO_->colspi()[0]; i++) {
                domains_[i]->print(out, i);
            }
        }

        if (auxiliary_domains_.size()) {
            fprintf(out, "  => Auxiliary Domains <=\n\n");
            for (int i = 0; i < L_AO_->colspi()[0]; i++) {
                auxiliary_domains_[i]->print(out, i);
            }
        }

        fflush(out);
    }
    double Local::er_metric()
    {
        int nso = L_AO_->rowspi()[0];
        int nmo = L_AO_->colspi()[0];

        boost::shared_ptr<Vector> iiii(new Vector("(ii|ii)", nmo));
        double* ip = iiii->pointer();

        boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
        boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
        const double* buffer = eri->buffer();
        double** Lp = L_AO_->pointer();

        for (int M = 0; M < basisset_->nshell(); M++) {
            for (int N = 0; N < basisset_->nshell(); N++) {
                for (int R = 0; R < basisset_->nshell(); R++) {
                    for (int S = 0; S < basisset_->nshell(); S++) {

                        int nM = basisset_->shell(M).nfunction();
                        int nN = basisset_->shell(N).nfunction();
                        int nR = basisset_->shell(R).nfunction();
                        int nS = basisset_->shell(S).nfunction();
                        int mstart = basisset_->shell(M).function_index();
                        int nstart = basisset_->shell(N).function_index();
                        int rstart = basisset_->shell(R).function_index();
                        int sstart = basisset_->shell(S).function_index();

                        eri->compute_shell(M,N,R,S);

                        for (int oM = 0, index = 0; oM < nM; oM++) {
                            for (int oN = 0; oN < nN; oN++) {
                                for (int oR = 0; oR < nR; oR++) {
                                    for (int oS = 0; oS < nS; oS++, index++) {
                                        for (int i = 0; i < nmo; i++) {
                                            ip[i] += Lp[oM + mstart][i] * Lp[oN + nstart][i] *
                                                buffer[index] *
                                                Lp[oR + rstart][i] * Lp[oS + sstart][i];
                                        }
                                    }}}}
                    }}}}

        double metric = 0.0;

        for (int i = 0; i < nmo; i++) {
            metric += ip[i];
        }

        return metric;
    }
    double Local::pm_metric()
    {
        int nso = L_AO_->rowspi()[0];
        int nmo = L_AO_->colspi()[0];
        int natom = basisset_->molecule()->natom();

        SharedMatrix Q = mulliken_charges(L_AO_);

        double metric = C_DDOT(nmo * (ULI) natom, Q->pointer()[0], 1, Q->pointer()[0], 1);

        return metric;
    }
    double Local::boys_metric()
    {
        double metric = 0.0;

        int nso = L_AO_->rowspi()[0];
        int nmo = L_AO_->colspi()[0];

        // Build dipole integrals
        boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
        boost::shared_ptr<OneBodyAOInt> Dint(integral->ao_dipole());

        std::vector<SharedMatrix > dipole;
        dipole.push_back(SharedMatrix(new Matrix("Dipole X", nso, nso)));
        dipole.push_back(SharedMatrix(new Matrix("Dipole Y", nso, nso)));
        dipole.push_back(SharedMatrix(new Matrix("Dipole Z", nso, nso)));

        Dint->compute(dipole);

        SharedMatrix XC(new Matrix("XC", nso, nmo));

        // X contribution
        C_DGEMM('N','N',nso,nmo,nso,1.0,dipole[0]->pointer()[0],nso,L_AO_->pointer()[0],nmo,0.0,XC->pointer()[0],nmo);

        for (int i = 0; i < nmo; i++) {
            double cont = C_DDOT(nso,&L_AO_->pointer()[0][i],nmo,&dipole[0]->pointer()[0][i],nmo);
            metric += cont * cont;
        }

        // Y contribution
        C_DGEMM('N','N',nso,nmo,nso,1.0,dipole[1]->pointer()[0],nso,L_AO_->pointer()[0],nmo,0.0,XC->pointer()[0],nmo);

        for (int i = 0; i < nmo; i++) {
            double cont = C_DDOT(nso,&L_AO_->pointer()[0][i],nmo,&dipole[1]->pointer()[0][i],nmo);
            metric += cont * cont;
        }

        // Z contribution
        C_DGEMM('N','N',nso,nmo,nso,1.0,dipole[2]->pointer()[0],nso,L_AO_->pointer()[0],nmo,0.0,XC->pointer()[0],nmo);

        for (int i = 0; i < nmo; i++) {
            double cont = C_DDOT(nso,&L_AO_->pointer()[0][i],nmo,&dipole[2]->pointer()[0][i],nmo);
            metric += cont * cont;
        }

        return metric;
    }
    SharedMatrix Local::C_USO()
    {
        return C_USO_;
    }
    SharedMatrix Local::C_AO()
    {
        if (!AO2USO_.get() || AO2USO_->nirrep() == 1)
            return C_USO_;

        int nao = AO2USO_->rowspi()[0];
        int nmo = 0;
        for (int h = 0; h < AO2USO_->nirrep(); h++)
            nmo += C_USO_->colspi()[h];

        SharedMatrix C = SharedMatrix(new Matrix("C (C1 Symmetry)",nao,nmo));
        double** Cp = C->pointer();

        int counter = 0;
        for (int h = 0; h < AO2USO_->nirrep(); h++) {
            int nsopi = AO2USO_->colspi()[h];
            int nmopi = C->colspi()[h];
            if (nsopi == 0 || nmopi == 0) continue;
            double** Ca = C_USO_->pointer(h);
            double** X = AO2USO_->pointer(h);

            C_DGEMM('N','N',nao,nmopi,nsopi,1.0,X[0],nsopi,Ca[0],nmopi,0.0,&Cp[0][counter],nmo);

            counter += nmopi;
        }
        return C;
    }
    SharedMatrix Local::L_AO()
    {
        return L_AO_;
    }
    SharedMatrix Local::AO2USO()
    {
        return AO2USO_;
    }
    void Local::localize(const std::string& algorithm, double conv)
    {
        if (boost::iequals(algorithm, "CHOLESKY"))
            localize_cholesky(conv);
        else if (boost::iequals(algorithm, "PM"))
            localize_pm(conv);
        else if (boost::iequals(algorithm, "BOYS"))
            localize_boys(conv);
        else if (boost::iequals(algorithm, "ER"))
            localize_er(conv);
        else
            throw PSIEXCEPTION("Localization algorithm not recognized");
    }
    void Local::localize_cholesky(double conv)
    {
        if (print_) {
            fprintf(outfile, "  ==> Localization: Cholesky <==\n\n");
        }

        L_AO_ = SharedMatrix(C_AO()->clone());
        L_AO_->set_name("L Cholesky (C1 Symmetry)");

        int nso = L_AO_->rowspi()[0];
        int nmo = L_AO_->colspi()[0];

        SharedMatrix D(new Matrix("D",nso,nso));
        double** Lp = L_AO_->pointer();
        double** Dp = D->pointer();

        C_DGEMM('N','T',nso,nso,nmo,1.0,Lp[0],nmo,Lp[0],nmo,0.0,Dp[0],nso);

        if (debug_)
            D->print();

        SharedMatrix L = D->partial_cholesky_factorize();

        if (debug_)
            L->print();

        int nmo2 = L->colspi()[0];
        double** Lp2 = L->pointer();

        if (nmo2 < nmo)
            throw PSIEXCEPTION("Local: Cholesky factor has smaller numerical rank than nmo!");

        for (int i = 0; i < nmo; i++) {
            C_DCOPY(nso, &Lp2[0][i], nmo2, &Lp[0][i], nmo);
        }

        if (debug_) {
            C_DGEMM('N','T',nso,nso,nmo,-1.0,Lp[0],nmo,Lp[0],nmo,1.0,Dp[0],nso);
            D->print(outfile, "Residual");
        }
    }

    void Local::localize_pm(double tightness)
    {
        L_AO_ = SharedMatrix(C_AO()->clone());
        L_AO_->set_name("L Pipek-Mezey (C1 Symmetry)");

        if (wavefunction_->nirrep() != 1) {
            fprintf(outfile, "Run in C1 Symmetry");
            abort();
        }

        int nocc = wavefunction_->doccpi()[0];

        fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);

        fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
        fprintf(outfile, "\t------------------------------------------------------------\n");

        Matrix U(nocc,nocc);
        Matrix V(nocc,nocc);
        Matrix VV(nocc,nocc);
        V.identity();

        int nfzc = wavefunction_->frzcpi()[0];
        int natom = wavefunction_->molecule()->natom();
        int nso = basisset_->nbf();
        double** C = C_USO_->pointer();
        double** S = S_->pointer();
        double alphalast = 1.0;

        std::vector <int> aostart(natom), aostop(natom);
        int offset = 0;
        aostart[0] = 0;
        for (int atom=0; atom < natom; atom++) {
            int nshell = basisset_->nshell_on_center(atom);
            for (int shell=0; shell < nshell; shell++)
                offset += basisset_->shell(basisset_->shell_on_center(atom, shell)).nfunction();
            aostop[atom] = offset;
        }

        for (int atom=0; atom < natom-1; atom++)
            aostart[atom + 1] = aostop[atom];

        for (int iter=0; iter < 100; iter++) {
            double P = 0.0;
            for(int i=nfzc; i < nocc; i++) {
                for(int A=0; A < natom; A++) {
                    double PiiA = 0.0;

                    for(int l=aostart[A]; l < aostop[A]; l++)
                        for(int k=0; k < nso; k++)
                            PiiA += C[k][i] * C[l][i] * S[k][l];

                    P += PiiA * PiiA;
                }
            }

            /* Compute 2x2 rotations for Pipek-Mezey localization */
            double alphamax = 0.0;
            for(int s=nfzc; s < nocc; s++) {
                for(int t=nfzc; t < s; t++) {

                    double Ast = 0.0;
                    double Bst = 0.0;

                    for(int A=0; A < natom; A++) {

                        double Pst = 0.0;
                        double Pss = 0.0;
                        double Ptt = 0.0;

                        for(int l=aostart[A]; l < aostop[A]; l++) {
                            for(int k=0; k < nso; k++) {
                                Pst += 0.5 * (C[k][s] * C[l][t] + C[l][s] * C[k][t]) * S[k][l];

                                Pss += C[k][s] * C[l][s] * S[k][l];

                                Ptt += C[k][t] * C[l][t] * S[k][l];
                            }
                        }

                        Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);
                        Bst += Pst * (Pss - Ptt);

                    } /* A-loop */

                    /* Compute the rotation angle */
                    double AB = Ast * Ast + Bst * Bst;
                    double alpha = 0.0;
                    if(fabs(AB) > 1.0e-12) {
                        double cos4a = -Ast/sqrt(AB);
                        alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
                    }
                    else alpha = 0.0;


                    /* Keep up with the maximum 2x2 rotation angle */
                    alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);

                    double Uss = cos(alpha);
                    double Utt = cos(alpha);
                    double Ust = sin(alpha);
                    double Uts = -Ust;

                    /* Now do the rotation */
                    for(int k=0; k < nso; k++) {
                        double Cks = C[k][s];
                        double Ckt = C[k][t];
                        C[k][s] = Uss * Cks + Ust * Ckt;
                        C[k][t] = Uts * Cks + Utt * Ckt;
                    }

                    U.identity();

                    U(s, s) = Uss;
                    U(t, t) = Utt;
                    U(s, t) = Ust;
                    U(t, s) = Uts;

                    VV.gemm(false, true, 1.0, V, U, 0.0);
                    V.copy(VV);
                } /* t-loop */
            } /* s-loop */

            double conv = fabs(alphamax) - fabs(alphalast);
            fprintf(outfile, "\t%4d  %20.10f  %20.10f  %6.3e\n", iter, P, alphamax, conv); fflush(outfile);

            if((iter > 2) && ((fabs(conv) < tightness) || alphamax == 0.0)) break;
            alphalast = alphamax;

        } /* iter-loop */

        /* Transform occupied orbital eigenvalues */
        SharedVector evals = wavefunction_->epsilon_a();
        double** F = block_matrix(nocc, nocc);
        for(int i=0; i < nocc; i++)
            for(int j=0; j < nocc; j++)
                for(int k=0; k < nocc; k++)
                    F[i][j] += V[0][k][i] * evals->get(k) * V[0][k][j];

        /* Compute a reordering array based on the diagonal elements of F */
        int* orb_order = init_int_array(nocc);
        int* orb_boolean = init_int_array(nocc);
        for(int i=0; i < nocc; i++) { orb_order[i] = 0;  orb_boolean[i] = 0; }

        int max = 0;
        for(int i=0,max=0; i < nocc; i++) /* First, find the overall maximum */
            if(fabs(F[i][i]) > fabs(F[max][max])) max = i;

        orb_order[0] = max;  orb_boolean[max] = 1;

        for(int i=1; i < nocc; i++) {
            max = 0;
            while(orb_boolean[max]) max++; /* Find an unused max */
            for(int j=0; j < nocc; j++)
                if((fabs(F[j][j]) >= fabs(F[max][max])) && !orb_boolean[j]) max = j;
            orb_order[i] = max; orb_boolean[max] = 1;
        }

        /* Now reorder the localized MO's according to F */
        double** Ctmp = block_matrix(nso,nocc);
        for(int i=0; i < nocc; i++)
            for(int j=0; j < nso; j++) Ctmp[j][i] = C[j][i];

        for(int i=0; i < nocc; i++) {
            int iold = orb_order[i];
            for(int j=0; j < nso; j++) C[j][i] = Ctmp[j][iold];
            evals->set(i, F[iold][iold]);
        }
        free_block(Ctmp);

//        fprintf(outfile, "\n\tPipek-Mezey Localized MO's (after sort):\n");
//        C_USO_->print();

        // Apparently transqt/transqt2 still use chkpt for obtaining the MO coefficients.
        // Dump the data to disk:
        double *values = evals->to_block_vector();
        double **vectors = C_USO_->to_block_matrix();
        chkpt_init(PSIO_OPEN_OLD);
        chkpt_wt_evals(values);
        chkpt_wt_scf(vectors);
        chkpt_close();
        delete[] vectors[0];
        delete[] vectors;
        delete[] values;
    }

    void Local::localize_boys(double conv)
    {
        throw FeatureNotImplemented("psi::Local","localize_boys",__FILE__,__LINE__);
        L_AO_ = SharedMatrix(C_AO()->clone());
        L_AO_->set_name("L Boys (C1 Symmetry)");

    }
    void Local::localize_er(double conv)
    {
        throw FeatureNotImplemented("psi::Local","localize_er",__FILE__,__LINE__);
        L_AO_ = SharedMatrix(C_AO()->clone());
        L_AO_->set_name("L Edmiston-Ruedenberg (C1 Symmetry)");

    }
    SharedMatrix Local::lowdin_charges(SharedMatrix C)
    {
        boost::shared_ptr<Molecule> molecule = basisset_->molecule();
        int nmo = C->colspi()[0];
        int nso = C->rowspi()[0];
        int natom = molecule->natom();
        SharedMatrix Q(new Matrix("Q: Gross Lowdin charges (nmo x natom)", nmo, natom));

        SharedMatrix XC(new Matrix("XC", nso, nmo));

        double** Cp  = C->pointer();
        double** Xp  = X_->pointer();
        double** XCp = XC->pointer();

        C_DGEMM('N','N',nso,nmo,nso,1.0,Xp[0],nso,Cp[0],nmo,0.0,XCp[0],nmo);

        double** Qp = Q->pointer();
        for (int i = 0; i < nmo; i++) {
            for (int m = 0; m < nso; m++) {
                int atom = basisset_->shell_to_center(basisset_->function_to_shell(m));
                Qp[i][atom] += XCp[m][i] * XCp[m][i];
            }
        }

        return Q;
    }
    SharedMatrix Local::mulliken_charges(SharedMatrix C)
    {
        boost::shared_ptr<Molecule> molecule = basisset_->molecule();
        int nmo = C->colspi()[0];
        int nso = C->rowspi()[0];
        int natom = molecule->natom();
        SharedMatrix Q(new Matrix("Q: Gross Mulliken charges (nmo x natom)", nmo, natom));

        SharedMatrix XC(new Matrix("XC", nso, nmo));

        double** Cp  = C->pointer();
        double** Sp  = S_->pointer();
        double** XCp = XC->pointer();

        C_DGEMM('N','N',nso,nmo,nso,1.0,Sp[0],nso,Cp[0],nmo,0.0,XCp[0],nmo);

        double** Qp = Q->pointer();
        for (int i = 0; i < nmo; i++) {
            for (int m = 0; m < nso; m++) {
                int atom = basisset_->shell_to_center(basisset_->function_to_shell(m));
                Qp[i][atom] += XCp[m][i] * Cp[m][i];
            }
        }

        return Q;
    }
    void Local::compute_boughton_pulay_domains(double Qcutoff)
    {
        boost::shared_ptr<Molecule> molecule = basisset_->molecule();
        int nmo = L_AO_->colspi()[0];
        int nso = L_AO_->rowspi()[0];
        int natom = molecule->natom();

        // Build Mulliken charges for current local coefficients
        Q_ = mulliken_charges(L_AO_);
        double** Qp = Q_->pointer();

        // Clear domains
        domains_.clear();
        auxiliary_domains_.clear();

        double** Sp = S_->pointer();
        double** Lp = L_AO_->pointer();

        // Premultiply SL
        SharedMatrix SL(new Matrix("SL", nso, nmo));
        double** SLp = SL->pointer();
        C_DGEMM('N','N',nso,nmo,nso,1.0,Sp[0],nso,Lp[0],nmo,0.0,SLp[0],nmo);

        // Build each domain
        for (int i = 0; i < nmo; i++) {
            std::set<int> total_atoms;

            // Rank atoms in this domain by Lowdin charges
            std::vector<std::pair<double, int> > charges;
            for (int A = 0; A < natom; A++) {
                charges.push_back(std::make_pair(Qp[i][A], A));
            }
            std::sort(charges.begin(), charges.end(), std::greater<std::pair<double, int> >());

            // Add atoms until domain becomes complete enough
            std::vector<int> atoms_in_domain;
            std::vector<int> funs_in_domain;
            for (int atom = 0; atom < natom; atom++) {
                // Add the current highest-charge atom to the domain
                int current_atom = charges[atom].second;
                atoms_in_domain.push_back(current_atom);

                // Add the functions from that atom to the funs_in_domain list
                for (int M = 0; M < basisset_->nshell(); M++) {
                    if (basisset_->shell_to_center(M) != current_atom) continue;
                    int nM = basisset_->shell(M).nfunction();
                    int mstart = basisset_->shell(M).function_index();
                    for (int om = 0; om < nM; om++) {
                        funs_in_domain.push_back(om + mstart);
                    }
                }

                // Figure out how many functions there are
                int nfun = funs_in_domain.size();

                // Temps
                SharedMatrix Smn(new Matrix("Smn", nfun, nfun));
                boost::shared_ptr<Vector> A(new Vector("A", nfun));
                double** Smnp = Smn->pointer();
                double* Ap = A->pointer();

                // Place the proper overlap elements
                for (int m = 0; m < nfun; m++) {
                    for (int n = 0; n < nfun; n++) {
                        Smnp[m][n] = Sp[funs_in_domain[m]][funs_in_domain[n]];
                    }
                }

                // Place the proper SL elements
                for (int m = 0; m < nfun; m++) {
                    Ap[m] = SLp[funs_in_domain[m]][i];
                }

                // Find the A vector
                int info1 = C_DPOTRF('L', nfun, Smnp[0], nfun);
                if (info1 != 0) {
                    throw PSIEXCEPTION("Local: Boughton Pulay Domains: Cholesky Failed!");
                }
                int info2 = C_DPOTRS('L', nfun, 1, Smnp[0], nfun, Ap, nfun);
                if (info2 != 0) {
                    throw PSIEXCEPTION("Local: Boughton Pulay Domains: Cholesky Solve Failed!");
                }

                // Compute the incompleteness metric
                double incompleteness = 1.0 - C_DDOT(nfun, Ap, 1, Ap, 1);
                if (incompleteness < 1.0 - Qcutoff)
                    break;
            }

            // rebuild a set<int>, which OrbtialDomain expects
            for (int A = 0; A < atoms_in_domain.size(); A++) {
                total_atoms.insert(atoms_in_domain[A]);
            }

            domains_.push_back(OrbitalDomain::buildOrbitalDomain(basisset_, total_atoms));
            if (auxiliary_.get() != NULL)
                auxiliary_domains_.push_back(OrbitalDomain::buildOrbitalDomain(auxiliary_, total_atoms));
        }
    }
    void Local::compute_polly_domains(double Qcutoff, double Rext, double Qcheck)
    {
        boost::shared_ptr<Molecule> molecule = basisset_->molecule();
        int nmo = L_AO_->colspi()[0];
        int nso = L_AO_->rowspi()[0];
        int natom = molecule->natom();

        // Build Lowdin charges for current local coefficients
        Q_ = lowdin_charges(L_AO_);
        double** Qp = Q_->pointer();

        // Clear domains
        domains_.clear();
        auxiliary_domains_.clear();

        // Build each domain
        for (int i = 0; i < nmo; i++) {
            std::set<int> charge_atoms;
            std::set<int> range_atoms;
            std::set<int> total_atoms;

            // Add atoms if they are greater than the charge cutoff
            for (int A = 0; A < natom; A++) {
                if (Qp[i][A] > Qcutoff) {
                    charge_atoms.insert(A);
                }
            }

            // Add atoms if they are inside the extended range parameter
            for (std::set<int>::const_iterator it = charge_atoms.begin();
                    it != charge_atoms.end(); it++) {
                Vector3 v = molecule->xyz(*it);
                for (int A = 0; A < natom; A++) {
                    if (v.distance(molecule->xyz(A)) < Rext / pc_bohr2angstroms) {
                        range_atoms.insert(A);
                    }
                }
            }

            // Merge the two lists
            for (std::set<int>::const_iterator it = charge_atoms.begin();
                    it != charge_atoms.end(); it++) {
                total_atoms.insert((*it));
            }
            for (std::set<int>::const_iterator it = range_atoms.begin();
                    it != range_atoms.end(); it++) {
                total_atoms.insert((*it));
            }

            // Check total charge to be sure the domain did not sneak by
            double charge_check = 0.0;
            for (std::set<int>::const_iterator it = total_atoms.begin();
                    it != total_atoms.end(); it++) {
                charge_check += Qp[i][(*it)];
            }

            // If the domain sneaks by, add all atoms to the domain
            if (charge_check < Qcheck) {
                total_atoms.clear();
                for (int A = 0; A < natom; A++) {
                    total_atoms.insert(A);
                }
            }

            domains_.push_back(OrbitalDomain::buildOrbitalDomain(basisset_, total_atoms));
            if (auxiliary_.get() != NULL)
                auxiliary_domains_.push_back(OrbitalDomain::buildOrbitalDomain(auxiliary_, total_atoms));
        }
    }
    OrbitalDomain::OrbitalDomain()
    {
    }
    OrbitalDomain::~OrbitalDomain()
    {
    }
    void OrbitalDomain::print(FILE* out, int label)
    {
        fprintf(out, "  => Orbital Domain %d: %ld Atoms <=\n\n", label, atoms_.size());

        fprintf(out, "    Atoms Local to Global:\n");
        for (int A = 0; A < atoms_local_to_global_.size(); A++) {
            fprintf(out, "    %4d -> %4d\n", A, atoms_local_to_global_[A]);
        }
        fprintf(out, "    \n");

        fprintf(out, "    Atoms Global to Local:\n");
        for (int A = 0; A < atoms_global_to_local_.size(); A++) {
            int a = atoms_local_to_global_[A];
            fprintf(out, "    %4d -> %4d\n", a, atoms_global_to_local_[a]);
        }
        fprintf(out, "    \n");

        fprintf(out, "    Shells Local to Global:\n");
        for (int A = 0; A < shells_local_to_global_.size(); A++) {
            fprintf(out, "    %4d -> %4d\n", A, shells_local_to_global_[A]);
        }
        fprintf(out, "    \n");

        fprintf(out, "    Shells Global to Local:\n");
        for (int A = 0; A < shells_global_to_local_.size(); A++) {
            int a = shells_local_to_global_[A];
            fprintf(out, "    %4d -> %4d\n", a, shells_global_to_local_[a]);
        }
        fprintf(out, "    \n");

        fprintf(out, "    Functions Local to Global:\n");
        for (int A = 0; A < functions_local_to_global_.size(); A++) {
            fprintf(out, "    %4d -> %4d\n", A, functions_local_to_global_[A]);
        }
        fprintf(out, "    \n");

        fprintf(out, "    Functions Global to Local:\n");
        for (int A = 0; A < functions_global_to_local_.size(); A++) {
            int a = functions_local_to_global_[A];
            fprintf(out, "    %4d -> %4d\n", a, functions_global_to_local_[a]);
        }
        fprintf(out, "    \n");
    }
    boost::shared_ptr<OrbitalDomain> OrbitalDomain::buildOrbitalDomain(boost::shared_ptr<BasisSet> basis, const std::set<int>& atoms)
    {
        boost::shared_ptr<OrbitalDomain> domain(new OrbitalDomain());
        domain->atoms_ = atoms;

        boost::shared_ptr<Molecule> molecule = basis->molecule();
        int natom = molecule->natom();
        int nso = basis->nbf();
        int nshell = basis->nshell();

        for (std::set<int>::const_iterator it = atoms.begin();
                it != atoms.end(); it++) {
            domain->atoms_local_to_global_.push_back((*it));
        }

        std::sort(domain->atoms_local_to_global_.begin(), domain->atoms_local_to_global_.end());

        for (int A = 0; A < atoms.size(); A++) {
            domain->atoms_global_to_local_[domain->atoms_local_to_global_[A]] = A;
        }

        int shell_counter = 0;
        int fun_counter = 0;
        for (int M = 0; M < nshell; M++) {
            int atom = basis->shell_to_center(M);
            if (!atoms.count(atom)) continue;

            domain->shells_local_to_global_.push_back(M);
            domain->shells_global_to_local_[M] = shell_counter++;

            int mstart = basis->shell(M).function_index();
            int nM = basis->shell(M).nfunction();

            for (int om = 0; om < nM; om++) {
                domain->functions_local_to_global_.push_back(om + mstart);
                domain->functions_global_to_local_[om + mstart] = fun_counter++;
            }
        }

        return domain;
    }

}} // Namespace psi

