using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using Rhino.PlugIns;
using System.Linq;
using Rhino.DocObjects;
using System.Drawing.Text;
using System.IO;

namespace MyPlugin.UtilityComps
{
    public class ImportPieces : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the StringConcatComponent class.
        /// </summary>
        public ImportPieces()
          : base("Import Pieces", "Import",
              "Concate two Strings",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("filePaths", "P", "First String to concatenate", GH_ParamAccess.list);
            pManager.AddPointParameter("Location Play", "L", "Import when clicked", GH_ParamAccess.item);
            //must be a toggle to avoid weird things
            pManager.AddBooleanParameter("import", "I", "Import when clicked", GH_ParamAccess.item);
            //set as a button to avoid hg weirdness
            pManager.AddBooleanParameter("reset", "r", "delete when done", GH_ParamAccess.item);
            
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("out", "O", "Length of concatanated String", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Geometry", "G", "Length of concatanated String", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        /// Importing Pieces from Junichiro Horikawa https://youtu.be/mljltzjgzyI?si=r9BN6aGQUVMlfNJv
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //go through the dir and import all smooth files give them the name of their indedx
            
            String stlPath = "";
            List<String> stlPaths = new List<String>();
            Point3d locPlay = new Point3d();
            Boolean import=false;
            Boolean reset=false;
            if (!DA.GetDataList(0,  stlPaths))return;
            if (!DA.GetData(1, ref locPlay)) return;
            if (!DA.GetData(2, ref import))return;
            if (!DA.GetData(3, ref reset)) return;

            //do work
            if (import)
            {
                
                    RhinoDoc.ActiveDoc.Objects.UnselectAll();
                foreach(string path in stlPaths) {
                    string cmd = "!_-Import " + "\"" + path + "\"" + " _Enter";
                    ///unselect all geometries
                    ///this line may be weird on Macs
                    //should be above for loop?
                    //import into rhino view via command line
                    DA.SetData(0, cmd);
                    Rhino.RhinoApp.RunScript(cmd, true);  //true shows the command in the command line
                    }
                    ///run command in rhino command which will actually make geometries in rhino 

                    /*bring geos into grasshopper*/
                    //get all geometry ids we have imported
                    var selectedObjs = RhinoDoc.ActiveDoc.Objects.GetSelectedObjects(false, false);
                    //access geometry property from selected object
                    geos = new List<GeometryBase>(); //geometryBase is versitile
                
                    foreach (RhinoObject selectedObj in selectedObjs)
                    {
                        var geo = selectedObj.Geometry;
                        var xf = Transform.Scale(Point3d.Origin+locPlay,1);
                        geos.Add(geo); //add to the list of geometries to display
                        geo.Transform(xf);
                        RhinoDoc.ActiveDoc.Objects.Delete(selectedObj, true); //delete the geometry in rhino but keeping it in gh
                }
            }

            ///remove all geometries from the list to delete them all from gh
            if (reset) {
                geos.Clear();
            }
            DA.SetDataList(1, geos);
            
            //clear the output list when reset is clicked
           
        }

        private List<GeometryBase> geos = new List<GeometryBase>(); // preserve geometies by storing them in a global var
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.import_icon;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("9449C992-222D-408F-84CD-3E1096C21A90"); }
        }
    }
}