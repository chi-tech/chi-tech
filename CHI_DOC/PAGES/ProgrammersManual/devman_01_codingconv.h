/**\page DevManCodeConventions Coding conventions

 \section sec0 Coding conventions

 In general we follow the
 <a href="https://google.github.io/styleguide/cppguide.html">Google C++ style guide</a>.
 We now will dictate a few specific items.
 \n\n
 Let us say it here explicitly, <B>all variable-, function-, method- names should be descriptive</B>,
 avoid abbreviation, rather sacrifice horizontal space for readability.

 \subsection devman1_sec0_1 File names

 - <B>Filenames</B> should be all lower case and may include "_" underscores. i.e.
   chi_tech_main.cc
 - Do not uses dashes `-` in file names
 - Do not use filenames that already exist in `/usr/include`
 - In general, make your filenames very specific. For example, use
   `http_server_logs.h` rather than `logs.h`. A very common case is to have a pair
   of files called, e.g., `foo_bar.h` and `foo_bar.cc`, defining a
   class called `FooBar`
 - Header files have a `.h` extension, c/c++ code files have a `.cc` extension.
 - Classes and Categorical functions should be contained in their own folder
   unless they are simple.
 - Folder depth shall be kept to a minimum and not exceed a depth of 3
   sub-levels unless very well justified, i.e. `ChiMath/Discretization/PWL`
 - Within a class folder the class declaration must be contained in its own
   header file. Exceptions to this are allowed for templates and in-line
   functions
 - When the member methods/functions become very large, consider functionally
   splitting them. I.e. use:
   - classname_00_constrdestr.cc for the constructor and destructor
   - classname_01_init.cc for an initialize call
   - classname_02_ops.cc for numerous small operations.

\image html "devman_Filenames.png" "Example 1" width=300px

 - If the class is really small then it should have a header file and an
   associated .cc file

  \image html "devman_Filenames1.png" "Example 2" width=350px

\subsection devman1_sec0_2 General code

 - <B>Variables</B> must use lower case and include "_" between words. i.e.
    num_faces.
 - <B>Namespaces</B> should be lower case and may include underscores
 - <B>FunctionNames</B> have no underscores and each word starts with a capital
   letter.
 - <B>Classes and structures</B> have each
   word start with a capital letter. For example chi_mesh::Vector is a structure
   called Vector and is defined in the namespace chi_mesh.

\code
int good_variable_name;    // Good variable name
int bvm;                   // Bad variable name
int DontUseThisName;       // Do not name variables with capital

void ThisIsAFunctionName(int parameter1)
{
  return;
}

class CoolDataStructure     // Standalone class declaration
{
public:
  int data;
};

namespace my_space            // Namespace definition
{
  class CoolerDataStructure;
};


class my_space::CoolerDataStructure
{
public:
  int data;
};

class my_space::CoolerDataStructure {
public:
  int data;
};
\endcode

\subsection devman1_sec0_3 Tabs, spaces and braces

In order to save horizontal space standard indents are to be 2 spaces instead
of 4 space equivalent tabs. Other than this the convention is flexible.\n\n

Curly braces, parentheses and block braces does not have a specific requirement
other than being used in a sense that is optimal with respect to readibility.\n\n

Generally we require that code span a maximum of 80 characters. This is not
a hard rule but will greatly enhance code reliability. Especially in
split-screen configurations.

\subsection devman1_sec0_4 Constants and Enumerators

Constants should either look like variable names with a "k" prepended to them
or look like macros. Constants can be defined in macros, enumerations or within
code segments.

\code
#define   PI       3.1415926535898
#define k_pi       3.1415926535898;

int DEFAULT_SETTING1   = 1;
int k_default_setting2 = 2;

enum UrlTableErrors {
  kOK = 0,
  kErrorOutOfMemory,
  kErrorMalformedInput,
};
enum AlternateUrlTableErrors {
  OK = 0,
  OUT_OF_MEMORY = 1,
  MALFORMED_INPUT = 2,
};
\endcode

\subsection devman1_sec0_5 Comments

 Be consistent in using comments. Comments should be clear and concise and
 convey the algorithm. Every class, structure or function should be
 supplied with doxygen comment styles at the top of the item.

\include "../../Modules/DiffusionSolver/lua/diffusion_execute.cc"

## Annotating Class-declerations:
<B>MAXIMIZE SCOPE MINIMIZE SPACE</B>.
Class declerations must always happen within header files and should aim to
provide the maximum amount of scope <B>within the minimum amount of vertical
space.</B> Member variables/structures should be on the top portion of the
class decleration and methods should be on the bottom.

\include "../../CHI_TECH/ChiMesh/SurfaceMesh/chi_surfacemesh.h"

 \section devman1_sec1 Header files

 - Avoid using header files for code definition. Use only for
   struct/class/function declaration
 - Use <B>#ifdef</B> guards on all header files
 - Avoid using forward declarations. Use header files instead.
 - Maintain a consistent order of include files.

 \subsection devman1_sec1_1 General rules
 Header files are primarily used to make function/class/struct declarations
 available to other compile units without the need to compile the associated
 definitions. Unless you have a good reason to do so, do not define code within
 header files. It results in extra compile times and is considered not good
 practice. Of course there are exceptions like templating but in general code
 definitions should be contained in .cc files, <B>header files are for code
 declaration only</B>.

 \subsection devman1_sec1_2 The #ifndef guard
 All header files should have #define guards to prevent multiple inclusion.
 The format of the symbol name should be _<CATEGORY>_<FILE>_H_.
 To guarantee uniqueness, they should be based on the full path in a project's
 source tree. For example, the file foo/src/bar/baz.h in project foo should
 have the following guard:

 \code
#ifndef FOO_BAR_BAZ_H_
#define FOO_BAR_BAZ_H_

...

#endif  // FOO_BAR_BAZ_H_
 \endcode

 \subsection devman1_sec1_3 Forward declarations
 Avoid using forward declarations where possible.
 Just #include the headers you need.


\section devman1_sec2 Namespace hell and typedefinitions

<B>Don't use the "using namespace" directive</B>. Rather rename
groups of namespaces. Use "typedef" where appropriate.


\code
// Forbidden -- This pollutes the namespace.
using namespace foo;

// Shorten access to some commonly used names in .cc files.
namespace baz = ::foo::bar::baz;

// Type definition for each cell there is a pair, number of faces and
// face areas for each face
typedef std::vector<std::pair<int, std::vector<double>>> CellFaceAreas;

CellFaceAreas cell_face_areas;
\endcode


\section devman1_sec3 Other general items

\subsection devman1_sec3_0 Line length

Each line of text in your code should be at most 80 characters long.
A line may exceed 80 characters if it is

- a comment line which is not feasible to split without harming
  readability, ease of cut and paste or auto-linking -- e.g.
  if a line contains an example command or a literal URL longer
  than 80 characters.
- a raw-string literal with content that exceeds 80 characters.
  Except for test code, such literals should appear near
  the top of a file.
- an include statement.
- a header guard


\subsection devman1_sec3_1 Boolean expressions

When you have a boolean expression that is longer than the
 standard line length, be consistent in how you break up the lines.

In this example, the logical AND operator is always at the end
 of the lines:

\code
if (this_one_thing > this_other_thing &&
    a_third_thing == a_fourth_thing &&
    yet_another && last_one) {
  ...
}
\endcode

\subsection devman1_sec3_2 Constructor Initializer lists

Constructor initializer lists can be all on one line or with
 subsequent lines indented four spaces.

The acceptable formats for initializer lists are:

\code
// When everything fits on one line:
MyClass::MyClass(int var) : some_var_(var) {
  DoSomething();
}

// If the signature and initializer list are not all on one line,
// you must wrap before the colon and indent 4 spaces:
MyClass::MyClass(int var)
    : some_var_(var), some_other_var_(var + 1) {
  DoSomething();
}

// When the list spans multiple lines, put each member on its own line
// and align them:
MyClass::MyClass(int var)
    : some_var_(var),             // 4 space indent
      some_other_var_(var + 1) {  // lined up
  DoSomething();
}

// As with any other code block, the close curly can be on the same
// line as the open curly, if it fits.
MyClass::MyClass(int var)
    : some_var_(var) {}
\endcode

\subsection devman1_sec3_3 Operators

Spacing between operators is flexible.

\code
// Assignment operators always have spaces around them.
x = 0;

// Other binary operators usually have spaces around them, but it's
// OK to remove spaces around factors.  Parentheses should have no
// internal padding.
v = w * x + y / z;
v = w*x + y/z;
v = w * (x + z);

// No spaces separating unary operators and their arguments.
x = -5;
++x;
if (x && !y)
  ...
\endcode

\subsection devman1_sec3_4 Templates and Casts

\code
// No spaces inside the angle brackets (< and >), before
// <, or between >( in a cast
std::vector<string> x;
y = static_cast<char*>(x);

// Spaces between type and pointer are OK, but be consistent.
std::vector<char *> x;
\endcode





 */