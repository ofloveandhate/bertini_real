using System.Reflection;
using System.Drawing;
using System.Collections.Generic;

namespace bertini_real
{
  
#pragma warning disable CA1416 // Validate platform compatibility
  public static class IconLoader
  {
    private static Dictionary<string, Bitmap> Cache { get; } = new Dictionary<string, Bitmap>();

    private readonly static int IconSize = 24;

    public static Bitmap GetIcon(string iconName)
    {
      try
      {
        if (Cache.TryGetValue(iconName, out var image))
            return image;

        var assembly = Assembly.GetExecutingAssembly();
        var names = assembly.GetManifestResourceNames(); // Enable for Debugging
        
        var imageStream = assembly.GetManifestResourceStream($"bertini_real.Properties.Resources.{iconName}");
        if (imageStream is null)
          return new Bitmap(IconSize, IconSize);
        
        Bitmap bitmap = new Bitmap(new Bitmap(imageStream), new Size(IconSize, IconSize));
        Cache.Add(iconName, bitmap);

        return bitmap;
      }
      catch
      {
        return new Bitmap(IconSize, IconSize);
      }
    }
  }
#pragma warning restore CA1416 // Validate platform compatibility
}