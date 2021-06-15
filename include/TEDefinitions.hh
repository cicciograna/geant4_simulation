#ifndef TEDEF_H
#define TEDEF_H 1

// World
static double worldLength    =   5000;
static double worldWidth     =   1000;

// Origin
static double originZpos     =   1120.;
static double originZoffset  =   0.;    
// static double originZoffset  = originZpos + vesselLength/2.;    // z-axis offset to place source in (0,0,0)

// Vessel
// static double vesselLength   = 2182.0;
static double vesselLength   = 2242.0;
static double vesselRadius   =  300.0;
static double vesselThick    =    5.0;
static double vesselZoffset  = vesselLength/2 + vesselThick + (originZpos - (vesselLength/2 + vesselThick));
// static double vesselThick    =   5.0;


// Target
static double targetLength   =    3.0;  // full length of MCP
static double targetRadius   =   60.0;  // Radius of MCP
static double targetZ        =  932. + targetLength;

// Electrodes
// First electrode is E000
static double eThicksmall    =    1.0; // 0.3;  // mm - Thickness of small electrodes
static double eThicklarge    =    1.0; // 0.6;  // mm - Thickness of large electrodes

static double eRiinit        =   66.0;  // mm - Inner radius of source region electrodes
static double eRoinit        =  120.0;  // mm - Outer radius of source region electrodes
static double eRoinitSpecial =   91.0;  // mm - Outer radius of innermost source region electrode

static double eRilarge       =  180.0;  // mm - Inner radius of large ion/ele region electrodes
static double eRolarge       =  240.0;  // mm - Outer radius of large ion/ele region electrodes

static double eRosmall       =  240.0;  // mm - Outer radius of small ion/ele region electrodes
static double eRismall       =  214.0;  // mm - Inner radius of small ion/ele region electrodes

static double eRifinal       =   60.0;  // mm - Outer radius of initial/final ion/ele region electrode

static double eCHheight      =   12.0;  // mm - Height of the coil holder
static double eCHRi          =   68.0;  // mm - Inner radius of the coil holder
static double eCHRo          =   91.0;  // mm - Outer radius of the coil holder


#endif
