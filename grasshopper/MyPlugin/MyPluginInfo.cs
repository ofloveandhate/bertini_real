﻿using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

//meta data for plugin

namespace MyPlugin
{
    public class MyPluginInfo : GH_AssemblyInfo
    {
        public override string Name => "MyPlugin";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("cf30f496-469e-4187-9e55-cfa50a0cd566");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";
    }
}