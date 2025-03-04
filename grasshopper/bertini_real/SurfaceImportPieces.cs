using System;
using System.Collections.Generic;
using Rhino;
using Rhino.DocObjects;

using Rhino.Geometry;
using Grasshopper.Kernel;


namespace bertini_real
{
    public class SurfaceImportStlPieces : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the StringConcatComponent class.
        /// </summary>
        public SurfaceImportStlPieces()
          : base("Surface Import STL Pieces", "ImportSurfSTLs",
              "brings in data to Grasshopper from stl files",
              "bertini_real", "Surface")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // the list of paths to load
            pManager.AddTextParameter("filePaths", "Ps", "The filepaths of the .stl files for the smooth pieces", GH_ParamAccess.list);

            //must be a toggle to avoid weird things
            pManager.AddBooleanParameter("import_toggle", "T", "This block runs when true", GH_ParamAccess.item);

            //set as a button to avoid hg weirdness
            pManager.AddBooleanParameter("reset_button", "B", "delete (and re-load if toggle true) when pressed", GH_ParamAccess.item);
            
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Commands", "Cs", "the commands that ran when importing", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Meshes", "Ms", "imported geometry from the stl files", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        /// Importing Pieces from Junichiro Horikawa https://youtu.be/mljltzjgzyI?si=r9BN6aGQUVMlfNJv
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //go through the dir and import all smooth files give them the name of their indedx
            
            List<String> stlPaths = new List<String>();
            List<String> commands = new List<String>();

            Boolean import_toggle=false;
            Boolean reset_button=false;

            if (!DA.GetDataList(0,  stlPaths))     return;
            if (!DA.GetData(1, ref import_toggle)) return;
            if (!DA.GetData(2, ref reset_button))  return;

            //do work
            if (import_toggle)
            {
                RhinoDoc.ActiveDoc.Objects.UnselectAll();
                foreach(string path in stlPaths) {
                    string cmd = "!_-Import " + "\"" + path + "\"" + " _Enter";
                    ///unselect all geometries
                    ///this line may be weird on Macs
                    //should be above for loop?
                    //import into rhino view via command line
                    
                    commands.Add(cmd);

                    Rhino.RhinoApp.RunScript(cmd, true);  //true shows the command in the command line
                    }
                    ///run command in rhino command which will actually make geometries in rhino 

                    /*bring geos into grasshopper*/
                    //get all geometry ids we have imported
                    var selectedObjs = RhinoDoc.ActiveDoc.Objects.GetSelectedObjects(false, false);
                    //access geometry property from selected object
                    geos = new List<GeometryBase>(); //geometryBase is versatile
                
                    foreach (RhinoObject selectedObj in selectedObjs)
                    {
                        var geo = selectedObj.Geometry;
                        var xf = Transform.Scale(Point3d.Origin,1);
                        geos.Add(geo); //add to the list of geometries to display
                        geo.Transform(xf);
                        RhinoDoc.ActiveDoc.Objects.Delete(selectedObj, true); //delete the geometry in rhino but keeping it in gh
                }
            }

            ///remove all geometries from the list to delete them all from gh
            if (reset_button) {
                geos.Clear();
            }
            DA.SetDataList(0, commands);
            DA.SetDataList(1, geos);
            
            //clear the output list when reset is clicked
           
        }

        private List<GeometryBase> geos = new List<GeometryBase>(); // preserve geometies by storing them in a global var
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("9449C992-222D-408F-84CD-3E1096C21A90"); }
        }
    }
}