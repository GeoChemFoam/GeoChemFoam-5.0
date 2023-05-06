/*---------------------------------------------------------------------------*\
License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_meshToMesh0_H
#define Foam_meshToMesh0_H

#include "fvMesh.H"
#include "HashTable.H"
#include "fvPatchMapper.H"
#include "scalarList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class indexedOctree;
class treeDataCell;

/*---------------------------------------------------------------------------*\
                         Class meshToMesh0 Declaration
\*---------------------------------------------------------------------------*/

class meshToMesh0
{
    // Private Data

        // Mesh References

        const fvMesh& fromMesh_;
        const fvMesh& toMesh_;

        //- fromMesh patch labels
        HashTable<label> fromMeshPatches_;

        //- toMesh patch labels
        HashTable<label> toMeshPatches_;

        //- Patch map
        HashTable<word> patchMap_;

        //- toMesh patch labels which cut the from-mesh
        HashTable<label> cuttingPatches_;

        //- Cell addressing
        labelList cellAddressing_;

        //- Boundary addressing
        labelListList boundaryAddressing_;

        //- Inverse-distance interpolation weights
        mutable unique_ptr<scalarListList> inverseDistanceWeightsPtr_;

        //- Inverse-volume interpolation weights
        mutable unique_ptr<scalarListList> inverseVolumeWeightsPtr_;

        //- Cell to cell overlap addressing
        mutable unique_ptr<labelListList> cellToCellAddressingPtr_;

        //- Overlap volume
        mutable scalar V_;


    // Private Member Functions

        //- Calculates mesh to mesh addressing pattern.
        //  For each cell from one mesh find the closest cell centre
        //  in the other mesh
        void calcAddressing();

        void cellAddresses
        (
            labelList& cells,
            const pointField& points,
            const fvMesh& fromMesh,
            const List<bool>& boundaryCell,
            const indexedOctree<treeDataCell>& oc
        ) const;

        void calculateInverseDistanceWeights() const;

        void calculateInverseVolumeWeights() const;

        void calculateCellToCellAddressing() const;

        const scalarListList& inverseDistanceWeights() const;

        const scalarListList& inverseVolumeWeights() const;

        const labelListList& cellToCellAddressing() const;


    // Private Static Data Members

        //- Direct hit tolerance
        static const scalar directHitTol;


public:

    // Declare name of the class and its debug switch
    ClassName("meshToMesh0");


    //- Enumeration specifying required accuracy
    enum order
    {
        MAP,
        INTERPOLATE,
        CELL_POINT_INTERPOLATE,
        CELL_VOLUME_WEIGHT
    };


    // Constructors

        //- Construct from the two meshes, the patch name map for the patches
        //  to be interpolated and the names of the toMesh-patches which
        //  cut the fromMesh
        meshToMesh0
        (
            const fvMesh& fromMesh,
            const fvMesh& toMesh,
            const HashTable<word>& patchMap,
            const wordList& cuttingPatchNames
        );

        //- Construct from the two meshes assuming there is an exact mapping
        //  between the patches
        meshToMesh0
        (
            const fvMesh& fromMesh,
            const fvMesh& toMesh
        );


    //- Destructor
    ~meshToMesh0() = default;


    //- Patch-field interpolation class
    class patchFieldInterpolator
    :
        public fvPatchFieldMapper
    {
        const labelList& directAddressing_;


    public:

        // Constructors

            //- Construct given addressing
            patchFieldInterpolator(const labelList& addr)
            :
                directAddressing_(addr)
            {}


        //- Destructor
        virtual ~patchFieldInterpolator() = default;


        // Member Functions

            label size() const
            {
                return directAddressing_.size();
            }

            bool direct() const
            {
                return true;
            }

            bool hasUnmapped() const
            {
                return false;
            }

            const labelList& directAddressing() const
            {
                return directAddressing_;
            }
    };


    // Member Functions

        // Access

            const fvMesh& fromMesh() const
            {
                return fromMesh_;
            }

            const fvMesh& toMesh() const
            {
                return toMesh_;
            }

            //- From toMesh cells to fromMesh cells
            const labelList& cellAddressing() const
            {
                return cellAddressing_;
            }

            //- Overlap volume
            scalar V() const
            {
                return V_;
            }


        // Interpolation

            //- Map field
            template<class Type, class CombineOp>
            void mapField
            (
                Field<Type>&,
                const Field<Type>&,
                const labelList& adr,
                const CombineOp& cop
            ) const;

            //- Interpolate field using inverse-distance weights
            template<class Type, class CombineOp>
            void interpolateField
            (
                Field<Type>&,
                const VolumeField<Type>&,
                const labelList& adr,
                const scalarListList& weights,
                const CombineOp& cop
            ) const;

            //- Interpolate field using inverse-volume weights
            template<class Type, class CombineOp>
            void interpolateField
            (
                Field<Type>&,
                const VolumeField<Type>&,
                const labelListList& adr,
                const scalarListList& weights,
                const CombineOp& cop
            ) const;


            //- Interpolate field using cell-point interpolation
            template<class Type, class CombineOp>
            void interpolateField
            (
                Field<Type>&,
                const VolumeField<Type>&,
                const labelList& adr,
                const vectorField& centres,
                const CombineOp& cop
            )const;


            //- Interpolate internal volume field
            template<class Type, class CombineOp>
            void interpolateInternalField
            (
                Field<Type>&,
                const VolumeField<Type>&,
                order=INTERPOLATE,
                const CombineOp& cop = eqOp<Type>()
            ) const;

            template<class Type, class CombineOp>
            void interpolateInternalField
            (
                Field<Type>&,
                const tmp<VolumeField<Type>>&,
                order=INTERPOLATE,
                const CombineOp& cop = eqOp<Type>()
            ) const;


            //- Interpolate volume field
            template<class Type, class CombineOp>
            void interpolate
            (
                VolumeField<Type>&,
                const VolumeField<Type>&,
                order=INTERPOLATE,
                const CombineOp& cop = eqOp<Type>()
            ) const;

            template<class Type, class CombineOp>
            void interpolate
            (
                VolumeField<Type>&,
                const tmp<VolumeField<Type>>&,
                order=INTERPOLATE,
                const CombineOp& cop = eqOp<Type>()
            ) const;


            //- Interpolate volume field
            template<class Type, class CombineOp>
            tmp<VolumeField<Type>> interpolate
            (
                const VolumeField<Type>&,
                order=INTERPOLATE,
                const CombineOp& cop = eqOp<Type>()
            ) const;

            template<class Type, class CombineOp>
            tmp<VolumeField<Type>> interpolate
            (
                const tmp<VolumeField<Type>>&,
                order=INTERPOLATE,
                const CombineOp& cop = eqOp<Type>()
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshToMesh0Templates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
