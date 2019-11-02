/*
 * This code is written by Tim Yuan
 * Created on 9th July, 2018
 * Last modified by Tim Yuan
 * Last modified on 9th July, 2018
 *
 * This code features density of z calculation
 * Use is at your own risk
 */

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */



#include <string>
#include <vector>
#include <iostream>
#include <gromacs/trajectoryanalysis.h>

using namespace gmx;

/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class AnalysisTemplate : public TrajectoryAnalysisModule
{
    public:
        AnalysisTemplate();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();
        double                           bin_width=0.01;
        double                           den_cutoff=16.5;
        double                            boxx,boxy,boxz;
        double                            ztmp, comz, boxxAvg, boxyAvg;
        int                              glob_nframes;
        int                              *count;
        int                              maxBin, minBin;
        int                              bin;
//        int         nr = sel.posCount();
//        real        frave = 0.0;




    private:
        class ModuleData;

        std::string                      fnDist_="z-den";
        double                           cutoff_;
        Selection                        refsel_;
        Selection                        sel_;
        AnalysisDataPlotSettings                  plotSettings_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
};


AnalysisTemplate::AnalysisTemplate()
    : cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


void
AnalysisTemplate::initOptions(IOptionsContainer          *options,
                              TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Calculate the density profile in z direction for a specified group",
        "(ie. OW in water)",
        "Take the z-COG of the second group specified as the z=0 plane\n\n\n",
        "Tim Yuan\n",
        "Sarupria Group\n",
        "Clemson University\n",
        "9th July, 2018\n"
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile()
                           .store(&fnDist_).defaultBasename("Iamsoooooopissed")
                           .description("Density profile in z-direction"));

    options->addOption(SelectionOption("basegroup")
                           .store(&refsel_).required()
                           .description("Reference group for the density profile"));

    options->addOption(SelectionOption("densitygroup")
                           .store(&sel_).required()
                           .description("Group to calculate the density "));


    options->addOption(DoubleOption("bin").store(&bin_width)
                           .description("Bin width for density calculation"));


    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void
AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         & /*top*/)
{
    nb_.setCutoff(cutoff_);

    data_.setColumnCount(0, 1);

    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);
    plotSettings_ = settings.plotSettings();

    maxBin = 0;
    minBin = 1000000;
    boxxAvg = 0;
    boxyAvg = 0;
    int size_count=1000/bin_width;
    count = new int[size_count]();

    std::cout<<"reference group picked: "<<refsel_.name()<<std::endl;
    std::cout<<"density group picked: "<<sel_.name()<<std::endl;
}

//---------------------------------------------------------------------------

/*
 * Define the parameters
 */

//---------------------------------------------------------------------------

void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{

    const Selection           &refsel = pdata->parallelSelection(refsel_);

    
    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, refsel);
    boxx=fr.box[0][0];
    boxy=fr.box[1][1];
    boxz=fr.box[2][2];
    comz = 0;   

    const Selection &sel   = pdata->parallelSelection(sel_);
    int              nr    = sel.posCount();
    int              nr_ref = refsel.posCount();

    //loop over reference group to determin the center of mass in the z direction

    for (int i = 0; i<nr_ref; ++i)
    {
        SelectionPosition p = refsel.position(i);
        ztmp = p.x()[2];
        ztmp = ztmp - boxz*floor(ztmp/boxz);
        comz += ztmp;
    }
    comz /= nr_ref;


    for (int i = 0; i < nr; ++i)
    {
        SelectionPosition p = sel.position(i);
        //without PBC 
        //ztmp = p.x()[2] - comz;
        //with PBC
        ztmp = (p.x()[2]-boxz*floor(p.x()[2]/boxz)) - comz;
        bin = (int)((ztmp+100.0)/bin_width);
        if (bin<=0)
        {
            std::cout<<"Oh!!No!!!!!Bin can not be less than zero"<<std::endl;
            exit(1);
        }
        else
        {
            count[bin]++;                
            if(bin>=maxBin)
            {
                maxBin = bin+1;
            }
            if(bin<minBin)
            {
                minBin = bin;
            }
        }
    }
    boxxAvg += boxx;
    boxyAvg += boxy;
}


void
AnalysisTemplate::finishAnalysis(int nframes)
{
    boxxAvg /= nframes;
    boxyAvg /= nframes;
    glob_nframes = nframes;
    std::cout<<"number of frames= "<<nframes<<std::endl;
    std::cout<<"Average box size in x = "<<boxxAvg<<std::endl;
    std::cout<<"Average box size in y = "<<boxyAvg<<std::endl;
    std::cout<<"bin width= "<<bin_width<<std::endl;
   
    AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(plotSettings_));
    plotm->setFileName(fnDist_);
    plotm->setTitle("Density Plot");
    plotm->setXAxisIsTime();
    plotm->setYLabel("Density in z-direction");
    plotm->appendLegend(sel_.name());
    plotm->appendLegend(refsel_.name());
    data_.addModule(plotm);
 
    std::string filename;
    filename =fnDist_.c_str();
    FILE * myfile;
    myfile = fopen(filename.c_str(), "a");

    for (int i=minBin; i<maxBin; ++i)
    {
        double out_bin=(i+0.5)*bin_width-100.0;

         fprintf(myfile, "%lf\t%lf\n",
             out_bin, (double)count[i]/((double)glob_nframes*(bin_width*boxxAvg*boxyAvg)));
    }
    fclose(myfile);
    delete [] count;
//delete the file and free the memory
}


void
AnalysisTemplate::writeOutput()
{
}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisTemplate>(argc, argv);
}
