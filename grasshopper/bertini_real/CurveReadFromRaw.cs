using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Windows.Forms;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace bertini_real
{
    public class CurveReadFromRaw : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CurveReadFromRaw class.
        /// </summary>
        public CurveReadFromRaw()
          : base("Curve Read From Raw", "CurveRead",
              "Read a curve from raw bertini_real output files",
              "bertini_real", "Curve")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            ///File to get data from. should be br_complete.json which is generated using write_piece() in python
            pManager.AddTextParameter("Directory", "D", "Directory containing a completed, sampled curve decomposition", GH_ParamAccess.item);
            Params.Input[0].Optional = false;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Vertices", "V", "the complete set of vertices of the curve.  Discards imaginary parts!!!", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Curve edges", "E", "the edges of the curve", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {  
            string folder = "";
            
            if (!DA.GetData(0, ref folder)) return;
            
            folder = Path.Combine(folder, "output_dim_1_comp_0");

            List<Point3d> vertices = read_vertices(folder);
            
            DA.SetDataList(0, vertices);
        }

        private List<Point3d> read_vertices(string folder) {



            List<Point3d> vertices = new List<Point3d>();
            var filenames = new List<string>();
            var filenames_lengths = new List<int>();
            var projections = new List<String>();

            using (var reader = new StreamReader(Path.Combine(folder, "V_samp.vertex")))
            {
                // Read first line and get number of vertices, projections, etc.
                var topline = reader.ReadLine().Split(' ');
                var numVertices = int.Parse(topline[0]);
                var numProjections = int.Parse(topline[1]);
                var numNaturalVarsInclHomCoord = int.Parse(topline[2]);
                var numVariables = numNaturalVarsInclHomCoord - 1;
                var numFilenames = int.Parse(topline[3]);

                reader.ReadLine(); // burn a line above the projection

                // Skip unused data
                for (int i = 0; i < (numProjections * numNaturalVarsInclHomCoord); i++)
                {
                    projections.Add(reader.ReadLine().Trim()); // burn a line above the projection
                }

                var line = reader.ReadLine(); // Skip a line

                // Get file names
                for (int i = 0; i < numFilenames; i++)
                {
                    filenames_lengths.Add(int.Parse(reader.ReadLine().Trim())); // Skip a line // this is the length of the filename.  ugh, remnants of C programming.
                    filenames.Add(reader.ReadLine().Trim());
                }

                // Read vertices
                for (int i = 0; i < numVertices; i++)
                {
                    line = reader.ReadLine();
                    while (string.IsNullOrWhiteSpace(line))
                    {
                        line = reader.ReadLine();
                    }
                    var numberOfVariables = int.Parse(line);

                    var temporaryPoint_complex = new List<Complex>();

                    for (int j = 0; j < numberOfVariables; j++)
                    {
                        var complexNum = reader.ReadLine().Split(' ');
                        var realPart = float.Parse(complexNum[0]);
                        var imaginaryPart = float.Parse(complexNum[1]);

                        temporaryPoint_complex.Add(new Complex(realPart, imaginaryPart));
                    }

                    var point = Dehomogenize(temporaryPoint_complex);

                    line = reader.ReadLine(); // read a line
                    var numProjectionValues = int.Parse(line);

                    var projection_values = new List<Complex>();
                    for (int l = 0; l < numProjectionValues; l++)
                    {
                        var complexNum = reader.ReadLine().Split(' ');
                        var realPart = float.Parse(complexNum[0]);
                        var imaginaryPart = float.Parse(complexNum[1]);
                        projection_values.Add(new Complex(realPart, imaginaryPart));
                    }

                    line = reader.ReadLine(); // read a line
                    var input_file_index = int.Parse(line.Trim());

                    line = reader.ReadLine(); // read a line
                    var vertexType = int.Parse(line.Trim());

                    var numPathNumbersEndingHere = int.Parse(reader.ReadLine().Trim());
                    var pathNumbersEndingHere = reader.ReadLine(); // i'm not going to bother to properly do things with this right now.  it's going to be a string of ws separated integers.

                    vertices.Add(ConvertToPoint3d(point));
                }
            }

            return vertices;
        }

        private Point3d ConvertToPoint3d(List<Complex> point){
            Point3d result = new Point3d();
            if (point.Count == 2){
                result.X = point[0].Real;
                result.Y = point[1].Real;
                result.Z = 0;
            }
            
           else if (point.Count == 3){
                result.X = point[0].Real;
                result.Y = point[1].Real;
                result.Z = point[2].Real;
            }

            else if (point.Count > 3){
                result.X = point[0].Real;
                result.Y = point[1].Real;
                result.Z = point[2].Real;
            }

            return result;
        }
        private List<Complex> Dehomogenize(List<Complex> point)
        {
            List<Complex> result = new List<Complex>();
            // Implement dehomogenize logic here

            for (int i = 1; i < point.Count; i++)
            {
                result.Add(point[i]/point[0]);
            }

            return result;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        /// 
        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("telephone.png");

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("233A833C-52BB-4D37-BADE-0C33B19B644D"); }
        }
    }
}