using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace BertiniReal
{
  public class BertiniRealInfo : GH_AssemblyInfo
  {
    public override string Name => "BertiniReal";

    //Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => null;

    //Return a short string describing the purpose of this GHA library.
    public override string Description => "";

    public override Guid Id => new Guid("2c394355-4a24-4e13-84a4-eeefac4b4054");

    //Return a string identifying you or your company.
    public override string AuthorName => "Caden Joergens and silviana amethyst";

    //Return a string representing your preferred contact details.
    public override string AuthorContact => "";
  }
}
