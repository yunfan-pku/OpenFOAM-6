/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ISATmanager.H"
#include "LUscalarMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class FuncType>
Foam::ISATmanager<FuncType>::ISATmanager(label in_n, label out_n, FuncType& func, const word& name_in,const dictionary& ISATDict)
:ISATDict_(ISATDict.subDict(name_in)),
tableTree_(in_n, out_n,readLabel(ISATDict_.lookup("maxNLeafs")),readLabel(ISATDict_.lookup("NtimeTag"))), 
pfunc(&func), 
epsilon_(1e-3), 
relepsilon_(0.0), 
scaleFactor_(out_n, out_n), 
scaleIn_(in_n, in_n), 
toleranceOut_(out_n, out_n), 
initToleranceIn_(in_n, in_n), 
init_elp_(in_n, in_n), 
timeSteps_(0),
checkInterval_(readLabel(ISATDict_.lookup("checkInterval"))),
maxDepthFactor_(readScalar(ISATDict_.lookup("maxDepthFactor"))),
nRetrieved_(0), 
nGrowth_(0), 
nAdd_(0), 
nCall_(0), 
treename_(name_in)
{
    for (int i = 0;i < out_n;i++)
        for (int j = 0;j < out_n;j++)
        {
            scaleFactor_[i][j] = 0;
            toleranceOut_[i][j] = 0;
        }
    for (int i = 0;i < in_n;i++)
        for (int j = 0;j < in_n;j++)
        {
            init_elp_[i][j] = 0;
            scaleIn_[i][j] = 0;
            initToleranceIn_[i][j] = 0;
        }
    for (int i = 0;i < out_n;i++)
    {
        scaleFactor_[i][i] = 1.0;
    }    
    for (int i = 0;i < in_n;i++)
    scaleIn_[i][i] = 1.0;
    scalarList toleranceOut_temp(ISATDict_.lookup("toleranceOut"));
    scalarList initToleranceIn_temp(ISATDict_.lookup("initToleranceIn"));
    scalarList scaleIn_temp(ISATDict_.lookup("scaleIn"));
    for (int i = 0;i < out_n;i++)
    {
        toleranceOut_[i][i] = 1.0 / toleranceOut_temp[i];
    }
    for (int i = 0;i < in_n;i++)
    {
        initToleranceIn_[i][i] = 1.0 / initToleranceIn_temp[i];
        scaleIn_[i][i] = scaleIn_temp[i];
    }
    /*
word Tname = "vaporfratree";
if (Tname == treename_)
{
    scaleFactor_[0][0] = 0;
}
*/
//scaleFactor_[0][0] = 1/10000.0;

//Info << "haha" << endl;
}

template<class FuncType>
Foam::ISATmanager<FuncType>::~ISATmanager()
{
    showPreformance();
    //tableTree_.balance(scaleIn_);
    //showPreformance();
    //tableTree_.balance();
    //word Tname = "vaporfratree";// "vaporfratree";
/*
    if (Tname == treename_)
        //if (tableTree_.size_ > 0)
    {
        scalarList in(3), out(1);
        in[0] = 0.699953;
        in[1] = 14980100;
        in[2] = 477.937;
        Info << "Tname=" << treename_ << endl;
        call(in, out);
        ISATleaf* pleaf;
        pleaf = search(in);
        Info << out[0] << ", point=" << pleaf->value() << ", data=" << pleaf->data() << endl;
        Info << "A= " << pleaf->A() << endl;
        Info << "EOA= " << pleaf->EOA() << endl;
        //pfunc->value(in, out);
        //Info << out[0] << endl;

    }
    */

}
template<class FuncType>
template <class ...Args>
void Foam::ISATmanager<FuncType>::add(const scalarList& value, Args... arg)
{
    //word Tname = "vaporfratree";

    ISATbinaryTree& T = tableTree_;
    ISATleaf* pleaf;
    if (T.size() == T.maxNLeafs())
    {
        //        Info<<"Deleting"<<endl;
        T.deleteLeaf(T.timeTagList().pop());
    }
    scalarList R(tableTree_.n_out());
    pfunc->value(value, R, arg...);
    pleaf = T.insertNewLeaf(value, R);
    pfunc->derive(value, pleaf->A(), arg...);
    scalarRectangularMatrix At(tableTree_.n_in_, tableTree_.n_in_);
    /*for (int i = 0;i < tableTree_.n_in_;i++)
        for (int j = 0;j < tableTree_.n_in_;j++)
        {
            At[i][j] = 100;
        }*/
        //Info<<pleaf->A()<<endl;

    //pleaf->EOA() = ((pleaf->A()) * scaleFactor_ * scaleFactor_ * (pleaf->A().T()) + init_elp_ * init_elp_) / (epsilon_ * epsilon_);//+ init_elp_
    pleaf->EOA() = scaleIn_ * initToleranceIn_ * initToleranceIn_ * scaleIn_ + scaleIn_ * (pleaf->A()) * toleranceOut_ * toleranceOut_ * (pleaf->A().T()) * scaleIn_;
    /*
    if (Tname == treename_)
    {
        Info << pleaf->EOA();
        if (pleaf->EOA()[0][0] < 1e-6)
        {
            FatalErrorInFunction << "Test:" <<pleaf->EOA()<< exit(FatalError);
        }
        Info << "Done"<<endl;
    }
    */
    // pleaf->EOA() = pleaf->EOA() + init_elp_ / (epsilon_ * epsilon_);
     //Info<<pleaf->EOA()<<endl;
    nAdd_++;
}
template<class FuncType>
Foam::ISATleaf* Foam::ISATmanager<FuncType>::search(const scalarList& value)
{
    ISATleaf* p;
    tableTree_.binaryTreeSearch(value, tableTree_.root_, p);
    return p;
}
template<class FuncType>
template <class ...Args>
void Foam::ISATmanager<FuncType>::call
(
    const Foam::scalarList& value, scalarList& out, Args... arg
)
{
    ISATleaf* pleaf;
    //word Tname = "psitree";
    //int aaa;
    //if (Tname == treename_)
    //    aaa = 5;
    if (!retrieve(value, out))
    {
        pleaf = search(value);

        if (pleaf == nullptr)
        {
            pfunc->value(value, out, arg...);
            add(value, arg...);
        }
        else
        {
            //if (Tname == treename_ && pleaf->value()[1] < 1.4982e+07 && pleaf->value()[1]>1.4980e+07 && pleaf->value()[2] < 550.48 && pleaf->value()[2]>550.478)
            //    aaa = 5;
            scalarList dvalue(value.size(), Zero);
            scalarList leafvalue(value.size(), Zero);
            for (int i = 0;i < tableTree_.n_in_;i++)
            {
                leafvalue[i] = pleaf->value()[i];
                dvalue[i] = value[i] - leafvalue[i];
            }
            for (int i = tableTree_.n_in_;i < value.size();i++)
            {
                leafvalue[i] = value[i];
            }
            pfunc->value(leafvalue + dvalue, out, arg...);
            scalarList out2(out.size());
            scalarList others = leafvalue - dvalue;
            bool flag_others = true;
            for (int i = 0;i < others.size();i++)
            {
                flag_others = flag_others && others[i] >= 0;
            }
            scalar sumother = 0;
            for (int i = 0;i < others.size() - 2;i++)
            {
                sumother += others[i];
            }
            flag_others = flag_others && sumother <= 1;
            if (!flag_others)
            {
                add(value, arg...);
            }
            else
            {
                pfunc->value(leafvalue - dvalue, out2, arg...);
                //Info<<out<<endl;

                //scalarList out2, value2;
                //value2 = 2 * pleaf->value_ - value;
                //value2[0] = value[0];
                //pfunc->value(value2, out2);
                if (!grow(pleaf, dvalue, out, out2))
                    add(value, arg...);
            }
        }
    }/*
    else {
        scalarList out2;
        pfunc->value(value, out2);
        if (fabs(out2[0] - out[0]) > 0.1)
        {
            pleaf = search(value);
            Info << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! I'm " << treename_ << ". The error is " << out[0] - out2[0] << endl;
        }
    }*/
    nCall_++;
    if (nCall_ % 1000 == 0)
    {
        showPreformance();
    }
}

template<class FuncType>
void Foam::ISATmanager<FuncType>::tablevalue
(
    const Foam::scalarList& value, scalarList& out
)
{
    ISATleaf* pleaf;
    pleaf = search(value);
    pleaf->eval(value, out);
}
template<class FuncType>
bool Foam::ISATmanager<FuncType>::retrieve
(
    const Foam::scalarList& value, scalarList& out
)
{
    bool retrieved(false);
    ISATleaf* plf;

    // If the tree is not empty
    if (tableTree_.size())
    {
        plf = search(value);

        // lastSearch keeps track of the chemPoint we obtain by the regular
        // binary tree search
        //lastSearch_ = phi0;
        if (plf->inEOA(value,scaleIn_))
        {
            retrieved = true;
        }
        // After a successful secondarySearch, phi0 store a pointer to the
        // found chemPoint
        /*
        else if (tableTree_.secondaryBTSearch(phiq, phi0))
        {
            retrieved = true;
        }
        else if (MRURetrieve_)
        {
           typename SLList
                <
                chemPointISAT<CompType, ThermoType>*
                >::iterator iter = MRUList_.begin();

            for (; iter != MRUList_.end(); ++iter)
            {
                phi0 = iter();
                if (phi0->inEOA(phiq))
                {
                    retrieved = true;
                    break;
                }
            }
        }
        */
    }
    // The tree is empty, retrieved is still false
    else
    {
        // There is no chempoints that we can try to grow
        //lastSearch_ = nullptr;
    }

    if (retrieved)
    {
        //phi0->increaseNumRetrieve();
        //scalar elapsedTimeSteps =
        //    this->chemistry_.timeSteps() - phi0->timeTag();

        // Raise a flag when the chemPoint has been used more than the allowed
        // number of time steps
        /*
        if (elapsedTimeSteps > chPMaxLifeTime_ && !phi0->toRemove())
        {
            cleaningRequired_ = true;
            phi0->toRemove() = true;
        }
        */
        /*
        lastSearch_->lastTimeUsed() = this->chemistry_.timeSteps();
        addToMRU(phi0);
        calcNewC(phi0, phiq, Rphiq);
        nRetrieved_++;
        */
        nRetrieved_++;
        tableTree_.eval(value, out);
        plf->increaseNumRetrieve();
        tableTree_.timeTagList().renew(plf->pTimeTagList());
        return true;
    }
    else
    {
        // This point is reached when every retrieve trials have failed
        // or if the tree is empty
        return false;
    }
}

template<class FuncType>
bool Foam::ISATmanager<FuncType>::grow
(
    ISATleaf* plf,
    const scalarList& dvalue,
    const scalarList& data1,
    const scalarList& data2
)
{
    // If the pointer to the chemPoint is nullptr, the function stops
    if (!plf)
    {
        return false;
    }

    // Raise a flag when the chemPoint used has been grown more than the
    // allowed number of time
    /*
    if (phi0->nGrowth() > maxGrowth_)
    {
        cleaningRequired_ = true;
        phi0->toRemove() = true;
        return false;
    }
    */

    // If the solution RphiQ is still within the tolerance we try to grow it
    // in some cases this might result in a failure and the grow function of
    // the chemPoint returns false
    scalarList ret1(data1.size()), ret2(data2.size());//, ret2(data.size());
    scalarList dvalue_m(plf->value().size());
    for (int i = 0;i < dvalue_m.size();i++)
    {
        dvalue_m[i] = dvalue[i];
    }
    //plf->eval(value2, ret2);
    plf->eval(plf->value() + dvalue_m, ret1);
    plf->eval(plf->value() - dvalue_m, ret2);
    //Info << distance(ret, data) << endl;
    /*
    bool rel_flag = true;
    for (int i = 0;i < tableTree_.n_out_;i++)
    {
        rel_flag = rel_flag && (fabs(ret[i] - data[i]) <= data[i] * relepsilon_);
    }
    */
    //if ((distance(ret, data) <= epsilon_ || distance(ret, data) <= relepsilon_ * norm(data)))
    //if (distance(ret1, data1) <= epsilon_ && distance(ret2, data2) <= epsilon_)
    if (normalized_distance(ret1, data1) <= 1.0 && normalized_distance(ret2, data2) <= 1.0)
    {
        plf->grow(plf->value() + dvalue_m,scaleIn_);
        tableTree_.timeTagList().renew(plf->pTimeTagList());
        nGrowth_++;
        return true;
    }
    // The actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
}
template<class FuncType>
double Foam::ISATmanager<FuncType>::distance(const scalarList& l, const scalarList& r)
{
    double sum = 0;
    for (int i = 0;i < l.size();i++)
        sum += sqr((l[i] - r[i]) * scaleFactor_[i][i]);
    return sqrt(sum);
}
template<class FuncType>
double Foam::ISATmanager<FuncType>::normalized_distance(const scalarList& l, const scalarList& r)
{
    double sum = 0;
    for (int i = 0;i < l.size();i++)
        sum += sqr((l[i] - r[i]) * toleranceOut_[i][i]);
    return sqrt(sum);
}
template<class FuncType>
double Foam::ISATmanager<FuncType>::norm(const scalarList& l)
{
    double sum = 0;
    for (int i = 0;i < l.size();i++)
        sum += sqr(l[i] * scaleFactor_[i][i]);
    return sqrt(sum);
}
template<class FuncType>
void Foam::ISATmanager<FuncType>::showPreformance() const
{
    Info << treename_ << ", ISAT performance: nCall=" << nCall_ << ", notCall=" << notCall << ", nRetrieved=" << nRetrieved_ << ", nGrowth=" << nGrowth_ << ", nAdd=" << nAdd_ << endl;
    Info <<"NtimeSteps:"<<timeSteps_<<",Treedepth:"<< tableTree_.depth()<<",Mindepth:"<< ceil(log2(tableTree_.size()))<<endl;

}

/*
template<class CompType, class ThermoType>
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::ISAT
(
    const dictionary& chemistryProperties,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryTabulationMethod<CompType, ThermoType>
    (
        chemistryProperties,
        chemistry
    ),
    tableTree_(chemistry, this->coeffsDict_),
    scaleFactor_(chemistry.nEqns() + ((this->variableTimeStep()) ? 1 : 0), 1),
    runTime_(chemistry.time()),
    chPMaxLifeTime_
    (
        this->coeffsDict_.lookupOrDefault("chPMaxLifeTime", INT_MAX)
    ),
    maxGrowth_(this->coeffsDict_.lookupOrDefault("maxGrowth", INT_MAX)),
    checkEntireTreeInterval_
    (
        this->coeffsDict_.lookupOrDefault("checkEntireTreeInterval", INT_MAX)
    ),
    maxDepthFactor_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "maxDepthFactor",
            (tableTree_.maxNLeafs() - 1)
           /(log(scalar(tableTree_.maxNLeafs()))/log(2.0))
        )
    ),
    minBalanceThreshold_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "minBalanceThreshold", 0.1*tableTree_.maxNLeafs()
        )
    ),
    MRURetrieve_(this->coeffsDict_.lookupOrDefault("MRURetrieve", false)),
    maxMRUSize_(this->coeffsDict_.lookupOrDefault("maxMRUSize", 0)),
    lastSearch_(nullptr),
    growPoints_(this->coeffsDict_.lookupOrDefault("growPoints", true)),
    nRetrieved_(0),
    nGrowth_(0),
    nAdd_(0),
    cleaningRequired_(false)
{
    if (this->active_)
    {
        dictionary scaleDict(this->coeffsDict_.subDict("scaleFactor"));
        label Ysize = this->chemistry_.Y().size();
        scalar otherScaleFactor = readScalar(scaleDict.lookup("otherSpecies"));
        for (label i=0; i<Ysize; i++)
        {
            if (!scaleDict.found(this->chemistry_.Y()[i].member()))
            {
                scaleFactor_[i] = otherScaleFactor;
            }
            else
            {
                scaleFactor_[i] =
                    readScalar
                    (
                        scaleDict.lookup(this->chemistry_.Y()[i].member())
                    );
            }
        }
        scaleFactor_[Ysize] = readScalar(scaleDict.lookup("Temperature"));
        scaleFactor_[Ysize + 1] = readScalar(scaleDict.lookup("Pressure"));
        if (this->variableTimeStep())
        {
            scaleFactor_[Ysize + 2] = readScalar(scaleDict.lookup("deltaT"));
        }
    }

    if (this->variableTimeStep())
    {
        nAdditionalEqns_ = 3;
    }
    else
    {
        nAdditionalEqns_ = 2;
    }

    if (this->log())
    {
        nRetrievedFile_ = chemistry.logFile("found_isat.out");
        nGrowthFile_ = chemistry.logFile("growth_isat.out");
        nAddFile_ = chemistry.logFile("add_isat.out");
        sizeFile_ = chemistry.logFile("size_isat.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::~ISAT()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::addToMRU
(
    chemPointISAT<CompType, ThermoType>* phi0
)
{
    if (maxMRUSize_ > 0 && MRURetrieve_)
    {
        // First search if the chemPoint is already in the list
        bool isInList = false;
        typename SLList <chemPointISAT<CompType, ThermoType>*>::iterator iter =
            MRUList_.begin();
        for ( ; iter != MRUList_.end(); ++iter)
        {
            if (iter() == phi0)
            {
                isInList = true;
                break;
            }
        }
        // If it is in the list, then move it to front
        if (isInList)
        {
            if (iter() != MRUList_.first())
            {
                // iter hold the position of the element to move
                MRUList_.remove(iter);

                // Insert the element in front of the list
                MRUList_.insert(phi0);
            }
        }
        else // chemPoint not yet in the list, iter is last
        {
            if (MRUList_.size() == maxMRUSize_)
            {
                MRUList_.remove(iter);
                MRUList_.insert(phi0);
            }
            else
            {
                MRUList_.insert(phi0);
            }
        }
    }
}


template<class CompType, class ThermoType>
void Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::calcNewC
(
    chemPointISAT<CompType, ThermoType>* phi0,
    const scalarField& phiq,
    scalarField& Rphiq
)
{
    label nEqns = this->chemistry_.nEqns(); // Species, T, p
    bool mechRedActive = this->chemistry_.mechRed()->active();
    Rphiq = phi0->Rphi();
    scalarField dphi(phiq-phi0->phi());
    const scalarSquareMatrix& gradientsMatrix = phi0->A();
    List<label>& completeToSimplified(phi0->completeToSimplifiedIndex());

    // Rphiq[i]=Rphi0[i]+A(i, j)dphi[j]
    // where Aij is dRi/dphi_j
    for (label i=0; i<nEqns-nAdditionalEqns_; i++)
    {
        if (mechRedActive)
        {
            label si = completeToSimplified[i];
            // The species is active
            if (si != -1)
            {
                for (label j=0; j<nEqns-2; j++)
                {
                    label sj = completeToSimplified[j];
                    if (sj != -1)
                    {
                        Rphiq[i] += gradientsMatrix(si, sj)*dphi[j];
                    }
                }
                Rphiq[i] +=
                    gradientsMatrix(si, phi0->nActiveSpecies())*dphi[nEqns - 2];
                Rphiq[i] +=
                    gradientsMatrix(si, phi0->nActiveSpecies() + 1)
                   *dphi[nEqns - 1];

                if (this->variableTimeStep())
                {
                    Rphiq[i] +=
                        gradientsMatrix(si, phi0->nActiveSpecies() + 2)
                       *dphi[nEqns];
                }

                // As we use an approximation of A, Rphiq should be checked for
                // negative values
                Rphiq[i] = max(0, Rphiq[i]);
            }
            // The species is not active A(i, j) = I(i, j)
            else
            {
                Rphiq[i] += dphi[i];
                Rphiq[i] = max(0, Rphiq[i]);
            }
        }
        else // Mechanism reduction is not active
        {
            for (label j=0; j<nEqns; j++)
            {
                Rphiq[i] += gradientsMatrix(i, j)*dphi[j];
            }
            // As we use a first order gradient matrix, Rphiq should be checked
            // for negative values
            Rphiq[i] = max(0, Rphiq[i]);
        }
    }
}


template<class CompType, class ThermoType>
bool Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::grow
(
    chemPointISAT<CompType, ThermoType>* phi0,
    const scalarField& phiq,
    const scalarField& Rphiq
)
{
    // If the pointer to the chemPoint is nullptr, the function stops
    if (!phi0)
    {
        return false;
    }

    // Raise a flag when the chemPoint used has been grown more than the
    // allowed number of time
    if (phi0->nGrowth() > maxGrowth_)
    {
        cleaningRequired_ = true;
        phi0->toRemove() = true;
        return false;
    }

    // If the solution RphiQ is still within the tolerance we try to grow it
    // in some cases this might result in a failure and the grow function of
    // the chemPoint returns false
    if (phi0->checkSolution(phiq, Rphiq))
    {
        return phi0->grow(phiq);
    }
    // The actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
}


template<class CompType, class ThermoType>
bool
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::cleanAndBalance()
{
    bool treeModified(false);

    // Check all chemPoints to see if we need to delete some of the chemPoints
    // according to the ellapsed time and number of growths
    chemPointISAT<CompType, ThermoType>* x = tableTree_.treeMin();
    while(x != nullptr)
    {
        chemPointISAT<CompType, ThermoType>* xtmp =
            tableTree_.treeSuccessor(x);

        scalar elapsedTimeSteps = this->chemistry_.timeSteps() - x->timeTag();

        if ((elapsedTimeSteps > chPMaxLifeTime_) || (x->nGrowth() > maxGrowth_))
        {
            tableTree_.deleteLeaf(x);
            treeModified = true;
        }
        x = xtmp;
    }

    MRUList_.clear();

    // Check if the tree should be balanced according to criterion:
    //  -the depth of the tree bigger than a*log2(size), log2(size) being the
    //      ideal depth (e.g. 4 leafs can be stored in a tree of depth 2)
    if
    (
        tableTree_.size() > minBalanceThreshold_
     && tableTree_.depth() >
        maxDepthFactor_*log(scalar(tableTree_.size()))/log(2.0)
    )
    {
        tableTree_.balance();
        treeModified = true;
    }

    // Return a bool to specify if the tree structure has been modified and is
    // now below the user specified limit (true if not full)
    return (treeModified && !tableTree_.isFull());
}


template<class CompType, class ThermoType>
void Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::computeA
(
    scalarSquareMatrix& A,
    const scalarField& Rphiq,
    const scalar rhoi,
    const scalar dt
)
{
    bool mechRedActive = this->chemistry_.mechRed()->active();
    label speciesNumber = this->chemistry_.nSpecie();
    scalarField Rcq(this->chemistry_.nEqns() + nAdditionalEqns_ - 2);
    for (label i=0; i<speciesNumber; i++)
    {
        label s2c = i;
        if (mechRedActive)
        {
            s2c = this->chemistry_.simplifiedToCompleteIndex()[i];
        }
        Rcq[i] = rhoi*Rphiq[s2c]/this->chemistry_.specieThermo()[s2c].W();
    }
    Rcq[speciesNumber] = Rphiq[Rphiq.size() - nAdditionalEqns_];
    Rcq[speciesNumber + 1] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 1];
    if (this->variableTimeStep())
    {
        Rcq[speciesNumber + 2] = Rphiq[Rphiq.size() - nAdditionalEqns_ + 2];
    }

    // Aaa is computed implicitly,
    // A is given by A = C(psi0, t0+dt), where C is obtained through solving
    // d/dt C(psi0, t) = J(psi(t))C(psi0, t)
    // If we solve it implicitly:
    // (C(psi0, t0+dt) - C(psi0, t0))/dt = J(psi(t0+dt))C(psi0, t0+dt)
    // The Jacobian is thus computed according to the mapping
    // C(psi0,t0+dt)*(I-dt*J(psi(t0+dt))) = C(psi0, t0)
    // A = C(psi0,t0)/(I-dt*J(psi(t0+dt)))
    // where C(psi0,t0) = I
    scalarField dcdt(speciesNumber + 2, Zero);
    this->chemistry_.jacobian(runTime_.value(), Rcq, dcdt, A);

    // The jacobian is computed according to the molar concentration
    // the following conversion allows the code to use A with mass fraction
    for (label i=0; i<speciesNumber; i++)
    {
        label si = i;

        if (mechRedActive)
        {
            si = this->chemistry_.simplifiedToCompleteIndex()[i];
        }

        for (label j=0; j<speciesNumber; j++)
        {
            label sj = j;
            if (mechRedActive)
            {
                sj = this->chemistry_.simplifiedToCompleteIndex()[j];
            }
            A(i, j) *=
              -dt*this->chemistry_.specieThermo()[si].W()
               /this->chemistry_.specieThermo()[sj].W();
        }

        A(i, i) += 1;
        // Columns for pressure and temperature
        A(i, speciesNumber) *=
            -dt*this->chemistry_.specieThermo()[si].W()/rhoi;
        A(i, speciesNumber + 1) *=
            -dt*this->chemistry_.specieThermo()[si].W()/rhoi;
    }

    // For the temperature and pressure lines, ddc(dTdt)
    // should be converted in ddY(dTdt)
    for (label i=0; i<speciesNumber; i++)
    {
        label si = i;
        if (mechRedActive)
        {
            si = this->chemistry_.simplifiedToCompleteIndex()[i];
        }

        A(speciesNumber, i) *=
            -dt*rhoi/this->chemistry_.specieThermo()[si].W();
        A(speciesNumber + 1, i) *=
            -dt*rhoi/this->chemistry_.specieThermo()[si].W();
    }

    A(speciesNumber, speciesNumber) = -dt*A(speciesNumber, speciesNumber) + 1;

    A(speciesNumber + 1, speciesNumber + 1) =
        -dt*A(speciesNumber + 1, speciesNumber + 1) + 1;

    if (this->variableTimeStep())
    {
        A(speciesNumber + 2, speciesNumber + 2) = 1;
    }

    // Inverse of (I-dt*J(psi(t0+dt)))
    LUscalarMatrix LUA(A);
    LUA.inv(A);

    // After inversion, lines of p and T are set to 0 except diagonal.  This
    // avoid skewness of the ellipsoid of accuracy and potential issues in the
    // binary tree.
    for (label i=0; i<speciesNumber; i++)
    {
        A(speciesNumber, i) = 0;
        A(speciesNumber + 1, i) = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::retrieve
(
    const Foam::scalarField& phiq,
    scalarField& Rphiq
)
{
    bool retrieved(false);
    chemPointISAT<CompType, ThermoType>* phi0;

    // If the tree is not empty
    if (tableTree_.size())
    {
        tableTree_.binaryTreeSearch(phiq, tableTree_.root(), phi0);

        // lastSearch keeps track of the chemPoint we obtain by the regular
        // binary tree search
        lastSearch_ = phi0;
        if (phi0->inEOA(phiq))
        {
            retrieved = true;
        }
        // After a successful secondarySearch, phi0 store a pointer to the
        // found chemPoint
        else if (tableTree_.secondaryBTSearch(phiq, phi0))
        {
            retrieved = true;
        }
        else if (MRURetrieve_)
        {
            typename SLList
            <
                chemPointISAT<CompType, ThermoType>*
            >::iterator iter = MRUList_.begin();

            for ( ; iter != MRUList_.end(); ++iter)
            {
                phi0 = iter();
                if (phi0->inEOA(phiq))
                {
                    retrieved = true;
                    break;
                }
            }
        }
    }
    // The tree is empty, retrieved is still false
    else
    {
        // There is no chempoints that we can try to grow
        lastSearch_ = nullptr;
    }

    if (retrieved)
    {
        phi0->increaseNumRetrieve();
        scalar elapsedTimeSteps =
            this->chemistry_.timeSteps() - phi0->timeTag();

        // Raise a flag when the chemPoint has been used more than the allowed
        // number of time steps
        if (elapsedTimeSteps > chPMaxLifeTime_ && !phi0->toRemove())
        {
            cleaningRequired_ = true;
            phi0->toRemove() = true;
        }
        lastSearch_->lastTimeUsed() = this->chemistry_.timeSteps();
        addToMRU(phi0);
        calcNewC(phi0, phiq, Rphiq);
        nRetrieved_++;
        return true;
    }
    else
    {
        // This point is reached when every retrieve trials have failed
        // or if the tree is empty
        return false;
    }
}


template<class CompType, class ThermoType>
Foam::label Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::add
(
    const scalarField& phiq,
    const scalarField& Rphiq,
    const scalar rho,
    const scalar deltaT
)
{
    label growthOrAddFlag = 1;
    // If lastSearch_ holds a valid pointer to a chemPoint AND the growPoints_
    // option is on, the code first tries to grow the point hold by lastSearch_
    if (lastSearch_ && growPoints_)
    {
        if (grow(lastSearch_, phiq, Rphiq))
        {
            nGrowth_++;
            growthOrAddFlag = 0;
            addToMRU(lastSearch_);
            //the structure of the tree is not modified, return false
            return growthOrAddFlag;
        }
    }

    // If the code reach this point, it is either because lastSearch_ is not
    // valid, OR because growPoints_ is not on, OR because the grow operation
    // has failed. In the three cases, a new point is added to the tree.
    if (tableTree().isFull())
    {
        // If cleanAndBalance operation do not result in a reduction of the tree
        // size, the last possibility is to delete completely the tree.
        // It can be partially rebuild with the MRU list if this is used.
        if (!cleanAndBalance())
        {
            DynamicList<chemPointISAT<CompType, ThermoType>*> tempList;
            if (maxMRUSize_>0)
            {
                // Create a copy of each chemPointISAT of the MRUList_ before
                // they are deleted
                typename SLList
                <
                    chemPointISAT<CompType, ThermoType>*
                >::iterator iter = MRUList_.begin();
                for ( ; iter != MRUList_.end(); ++iter)
                {
                    tempList.append
                    (
                        new chemPointISAT<CompType, ThermoType>(*iter())
                    );
                }
            }
            tableTree().clear();

            // Pointers to chemPoint are not valid anymore, clear the list
            MRUList_.clear();

            // Construct the tree without giving a reference to attach to it
            // since the structure has been completely discarded
            chemPointISAT<CompType, ThermoType>* nulPhi = 0;
            forAll(tempList, i)
            {
                tableTree().insertNewLeaf
                (
                    tempList[i]->phi(),
                    tempList[i]->Rphi(),
                    tempList[i]->A(),
                    scaleFactor(),
                    this->tolerance(),
                    scaleFactor_.size(),
                    nulPhi
                );
                deleteDemandDrivenData(tempList[i]);
            }
        }

        // The structure has been changed, it will force the binary tree to
        // perform a new search and find the most appropriate point still stored
        lastSearch_ = nullptr;
    }

    // Compute the A matrix needed to store the chemPoint.
    label ASize = this->chemistry_.nEqns() + nAdditionalEqns_ - 2;
    scalarSquareMatrix A(ASize, Zero);
    computeA(A, Rphiq, rho, deltaT);

    tableTree().insertNewLeaf
    (
        phiq,
        Rphiq,
        A,
        scaleFactor(),
        this->tolerance(),
        scaleFactor_.size(),
        lastSearch_ // lastSearch_ may be nullptr (handled by binaryTree)
    );
    if (lastSearch_ != nullptr)
    {
        addToMRU(lastSearch_);
    }
    nAdd_++;

    return growthOrAddFlag;
}


template<class CompType, class ThermoType>
void
Foam::chemistryTabulationMethods::ISAT<CompType, ThermoType>::writePerformance()
{
    if (this->log())
    {
        nRetrievedFile_()
            << runTime_.timeOutputValue() << "    " << nRetrieved_ << endl;
        nRetrieved_ = 0;

        nGrowthFile_()
            << runTime_.timeOutputValue() << "    " << nGrowth_ << endl;
        nGrowth_ = 0;

        nAddFile_()
            << runTime_.timeOutputValue() << "    " << nAdd_ << endl;
        nAdd_ = 0;

        sizeFile_()
            << runTime_.timeOutputValue() << "    " << this->size() << endl;
    }
}
*/

// ************************************************************************* //
