using Rhino.Collections;
using Rhino.Geometry;
using System;


/*
 * A Sanity File to hold all the classes used in the code. 
 * 
 * Data used in transformConnectors to store data from the JSON
 * PieceData hold data on individual piece objects which has been parsed in TransformConnectors
 * The rest are used to create the parts of a plug in Connectors.cs or PositivePlugComponent.cs
 * plugBody create the main section of the plug with the tapered top
 * plugtabs create the cutout box or wedges for the plug
 * 
 * 
 */
namespace BertiniReal.UtilityComps
{
    /// <summary>
    /// Organize and store data from br_surf_piece_data.json file
    /// NOTE: The JSON no longer has piece_indicies, but now piece_names which are file name strings 
    /// </summary>
    /// <see cref="TransformConnectors.cs"/>
    public class Data
    {
        public int[][] piece_indices { get; set; } //this will need to change
        public int[][] singularities_on_pieces { get; set; }
        public double[][] sing_directions { get; set; }
        public double[][] sing_locations { get; set; }
        public int[][] parities { get; set; }
    }

    /// <summary>
    /// Store data parsed from the Data class by peice
    /// </summary>
    /// <see cref="TransformConnectors.cs"/>
    public class PieceData
    {
        public int piece_index { get; set; }
        public int[] indices { get; set; }
        public int[] singsOnPiece { get; set; }
        public Vector3d[] directions { get; set; }
        public Vector3d[] locations { get; set; }
        public int[] parities { get; set; }

        public PieceData()
        {

        }
    }

    /* Used to Create the different parts of a positive plug */
    public class PlugBody
    {
        public double taperedR { get; set; }
        public double plugTaperLength { get; set; }
        public double plugR { get; set; }
        public double bodyHeight { get; set; }
        public PlugBody(double taperedR, double plugTaperLength, double plugR, double bodyHeight)
        {
            this.taperedR = taperedR;
            this.plugTaperLength = plugTaperLength;
            this.plugR = plugR;
            this.bodyHeight = bodyHeight;
        }
        public Brep plugBodyGeo()
        {
            Point3d p0 = new Point3d(0, 0, 0);
            Point3d p1 = new Point3d(taperedR, 0, 0); //origin to taper radius
            Point3d p2 = new Point3d(plugR, 0, plugTaperLength); //transition from taper to body
            Point3d p3 = new Point3d(plugR, 0, plugTaperLength + bodyHeight); //to top of plug
            Point3d p4 = new Point3d(0, 0, plugTaperLength + bodyHeight); //top on axis

            Polyline polyline = new Polyline(new Point3dList(p0, p1, p2, p3));
            Line axis = new Line(p0, p4);
            Brep geo = RevSurface.Create(polyline, axis).ToBrep();
            return geo.CapPlanarHoles(0.001);
        }
    }
    /** A class to create all the piceses needed for the tabs on a positive plug*/
    public class PlugTabs
    {
        public double wedgeHeight { get; set; }
        public double wedgeLength { get; set; }
        public double wedgeWidth { get; set; }
        public double tabCutoutThickness { get; set; }
        public double tabCutoutDepth { get; set; }
        public double tabThickness { get; set; }

        public PlugTabs(double wedgeHeight, double wedgeLength, double wedgeWidth, double tabCutoutThickness, double tabCutoutDepth, double tabThickness)
        {
            this.wedgeHeight = wedgeHeight;
            this.wedgeLength = wedgeLength;
            this.wedgeWidth = wedgeWidth;
            this.tabCutoutThickness = tabCutoutThickness;
            this.tabCutoutDepth = tabCutoutDepth;
            this.tabThickness = tabThickness;
        }

        /** Create a box to subtract from the plug body
         * Will need to mirror to subtract from both sides
         */
        public Brep cutoutBox(double plugR, Plane b, double eps)
        {
            ///Creaft the box of the profile
            double ell = plugR - tabThickness - tabCutoutThickness / 2;
            Interval intervalX = new Interval(-plugR, plugR + 1);
            Interval intervalY = new Interval(-tabCutoutThickness / 2, tabCutoutThickness / 2);
            Interval intervalZ = new Interval(-(1 + eps + tabCutoutDepth) / 2, (eps + tabCutoutDepth) / 2);

            ///create the physical box from the bounds
            Brep cutoutBox = new Box(b, intervalX, intervalY, intervalZ).ToBrep();
            ///Transform to the center to easy transformation
            var xf = Transform.Translation(0, ell, tabCutoutDepth / 2);
            cutoutBox.Transform(xf);

            return cutoutBox;
        }

        /** Create a wedge which goes on the tab of the plug. Will need to be mirrored */
        public Brep rawWedge()
        {
            ///Make a triangle for wedge profile with right angle at the origin using a series of points
            //height(along zaxis) = plug_wedge_h, the height of the wedge
            //length(along yaxis) = plug_wedge_l+plug_tab_thickness, the lenght increases as the gap nears the center
            Point3d p0 = new Point3d(0, 0, 0);
            Point3d p1 = new Point3d(0, 0, -wedgeHeight);
            Point3d p2 = new Point3d(0, wedgeLength + tabThickness, 0);
            Vector3d extrudePath = new Vector3d(wedgeWidth, 0, 0); //path to extrude along = width of wedge

            //build the wedge profile in the 4th quad by connecting the points 
            PolylineCurve profile = new PolylineCurve(new Point3dList(p1, p0, p2, p1));

            ///Create the wedge by extruding the profile to plug_wedge_w
            Brep wedge = Extrusion.Create(profile, wedgeWidth, true).ToBrep();

            ///and then re-center the width on the x-axis to allow for easy transformation
            var xf = Transform.Translation(wedgeWidth / 2, 0, 0);
            wedge.Transform(xf);

            return wedge;
        }
    }

}
