import subprocess
import os
import sys
import textwrap

print("Creating InputParamsDocs")
main = open("doc/generated_files/test.txt", "r")

obj_list = []
object_dict = {}
object_dict["optional_params"] = []
object_dict["required_params"] = []
param_dict = {}

def EverythingAfter(line_string : str, word : str):
    start = line_string.find(word)
    end = start + len(word)

    return line_string[end:]

lines = main.readlines()
for line in lines:
    words = line.split()

    if len(words) <= 1:
        continue

    if words[1] == "OBJECT_BEGIN":
        object_dict["name"] = words[2]

    if words[1] == "OBJECT_END":
        obj_list.append(object_dict)
        object_dict = {}
        object_dict["optional_params"] = []
        object_dict["required_params"] = []

    if words[1] == "PARAM_BEGIN":
        param_dict["name"] = words[2]

    if words[1] == "PARAM_END":
        if param_dict["TAG"] == "OPTIONAL":
            object_dict["optional_params"].append(param_dict)
        else:
            object_dict["required_params"].append(param_dict)
        param_dict = {}

    if words[1] == "TYPE":
        param_dict[words[1]] = words[2]

    if words[1] == "TAG":
        param_dict[words[1]] = words[2]

    if words[1] == "DOC_STRING":
        param_dict[words[1]] = EverythingAfter(line, words[1]).strip()

    if words[1] == "DEFAULT_VALUE":
        param_dict[words[1]] = EverythingAfter(line, words[1]).strip()

    if words[1] == "CONSTRAINTS":
        param_dict[words[1]] = EverythingAfter(line, words[1]).strip()

main.close()

if not os.path.exists("doc/generated_files"):
    os.mkdir("doc/generated_files")

if not os.path.exists("doc/generated_files/input_docs"):
    os.mkdir("doc/generated_files/input_docs")

dropdn_list_script1 = '''
{
    var coll = document.getElementsByClassName("droppy1");
    var i;
    
    for (i = 0; i < coll.length; i++) {
      coll[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var content = this.parentElement.nextElementSibling;
        if (content.style.display === "block") {
          content.style.display = "none";
        } else {
          content.style.display = "block";
        }
      });
    }
}'''

dropdn_list_script2 = '''
{
    var coll = document.getElementsByClassName("droppy2");
    var i;
    
    for (i = 0; i < coll.length; i++) {
      coll[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var content = this.nextElementSibling;
        if (content.style.display === "block") {
          content.style.display = "none";
          this.childNodes[0].innerText = "►"
        } else {
          content.style.display = "block";
          this.children[0].innerText = "▼"
        }
      });
    }
}'''

#background-color: #313335;
#width: 100%;
# padding: 2px 2px;

button_style = '''
"
background-color: #edf0f5;
color: #9373A5;
cursor: pointer;
border: none;
text-align: left;
outline: none;
font-size: 15px;
margin: 2px 2px;
border: 2px solid #687372;
display: block;
width: 100%;
vertical-align: middle;
padding-top: 4px;
padding-bottom: 0px;
"
'''
#<span class="arrow" style="padding-left: 0px;">▼
#<span class="arrow" style="padding-left: 16px;">►</span>

''' This function writes HTML formatted parameters with
    dropdown boxes.'''
def WriteHTMLParameters(file, params):
    file.write('\htmlonly\n\n')
    file.write('<div style="display: block;">\n')
    droppy = "droppy2"
    for param in params:
        file.write('<button type="button" class="'+droppy+'" ' + \
                   'style=' + button_style +'>' + \
                   '<span class="arrow" style="padding-left: 0px;">►' + \
                   '</span><TT><B>' + \
                   param["name"] + '</B></TT></button>\n')
        # if droppy == "droppy1":
        #     droppy = "droppy2"
        file.write('<div class="content" style="display: none;">\n')
        file.write('  <p><I>type=</I><span style="color: blue;"><TT>' + \
                   param["TYPE"] + '</TT></span>. ')
        file.write(param["DOC_STRING"] + "</p>\n")
        if "DEFAULT_VALUE" in param:
            file.write('  <p>Default value: <TT style="color:grey">' + \
                       param["DEFAULT_VALUE"] + \
                       "</TT></p>\n")
        if "CONSTRAINTS" in param:
            file.write("  <p>Allowable values: `" + param["CONSTRAINTS"] + \
                       "`</p>\n")
        file.write('</div>\n\n')

    file.write('</div>\n')
    file.write('\endhtmlonly\n\n')

for obj in obj_list:
    reduced_name = obj["name"].replace(":", "_")
    file = open("doc/generated_files/input_docs/" + reduced_name+".h", "w")

    file.write("/** \\addtogroup " + reduced_name + "\n\n\n")

    file.write(textwrap.dedent('''
    ## Example usage:
    Create this object:
    \code
    params = 
    {
        param_name1 = value1,
        param_name2 = value2,
        --etc.
    }\n'''))
    file.write(obj["name"].replace("::",".")+".Create(params)\n")

    file.write("\endcode\n")

    if len(obj["required_params"]) > 0:
        file.write("## Required Input parameters\n")
        WriteHTMLParameters(file, obj["required_params"])

    if len(obj["optional_params"]) > 0:
        file.write("## Optional Input parameters\n")
        WriteHTMLParameters(file, obj["optional_params"])

    file.write('\htmlonly\n\n')
    file.write('<script>\n')
    file.write(dropdn_list_script1)
    file.write(dropdn_list_script2)
    file.write('</script>\n\n')
    file.write('\endhtmlonly\n\n')

    file.write("*/")
    file.close()

print("Done creating InputParamsDocs")