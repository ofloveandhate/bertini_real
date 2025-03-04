using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace bertini_real
{
  public class bertini_realInfo : GH_AssemblyInfo
  {
    public override string Name => "bertini_real Info";

    //Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => null;

    //Return a short string describing the purpose of this GHA library.
    public override string Description => "";

    public override Guid Id => new Guid("09237482-2356-4177-8520-94ace80b49d3");

    //Return a string identifying you or your company.
    public override string AuthorName => "";

    //Return a string representing your preferred contact details.
    public override string AuthorContact => "";

    //Return a string representing the version.  This returns the same version as the assembly.
    public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
  }
}