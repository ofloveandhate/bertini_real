using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace bertini_real
{
  public class bertini_realInfo : GH_AssemblyInfo
  {
    public override string Name => "bertini_real";

    //Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => null;

    //Return a short string describing the purpose of this GHA library.
    public override string Description => "Work with real algebraic curves and surfaces computed with bertini_real";

    public override Guid Id => new Guid("09237482-2356-4177-8520-94ace80b49d3");

    //Return a string identifying you or your company.
    public override string AuthorName => "silviana amethyst and her many students and collaborators";

    //Return a string representing your preferred contact details.
    public override string AuthorContact => "amethyst@mpi-cbg.de";

    //Return a string representing the version.  This returns the same version as the assembly.
    public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
  }
}