#ifndef GNUPLOT_H
#define GNUPLOT_H

#include <cassert>
#include <fstream>
#include <string>

/*! \class GnuPlot
 * \brief Plot output data with Gnuplot.
 *
 * This class constructs a command line and passes it to gnuplot via \c std::system.
 */
class GnuPlot
{
public:
  GnuPlot(std::string _filename, std::ofstream &_output_file) :
    filename(_filename), output_file(_output_file)
  {
    output_file.open(filename);

    if (!output_file.is_open())
      throw std::runtime_error("file could not be opened");
  }

  ~GnuPlot()
  {
    output_file.close();
  }

  void plot_with_lines(size_t dim, std::string style = "lines",
                       bool plot_3d = false)
  {
    std::string cmd;

    switch (dim)
      {
      case 1:
        cmd = "gnuplot -p -e \"plot '"
            + filename
            + "' using 1:2 with "
            + style + "\"";
        break;

      case 2:
        if (plot_3d)
          cmd = "gnuplot -p -e \"splot '"
              + filename
              + "' using 1:2:3 with "
              + style + "\"";
        else
          cmd = "gnuplot -p -e \"plot '"
              + filename
              + "' using 1:2 with "
              + style + ", '"
              + filename
              + "'using 1:3 with "
              + style + "\"";
        break;

      default:
        throw std::invalid_argument("dimension >2 not implemented");
      }

    std::system(cmd.c_str());
  }

private:
  std::string filename;
  std::ofstream &output_file;
};

#endif // GNUPLOT_H
