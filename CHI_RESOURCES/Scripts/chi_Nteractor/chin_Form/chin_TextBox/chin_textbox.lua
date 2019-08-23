

TextBoxClass={}
TextBoxClass.__index = TextBoxClass
TextBoxCount = 0;

dofile(chinTextBoxDir.."chin_textbox_00_constrdestr.lua")
dofile(chinTextBoxDir.."chin_textbox_01_functions.lua")


--#################################
--         Design Notes
--#################################
--Textboxes consist of:
--Outline
--Back pane (with material)
--Label
--Input Label
--Cursor line
--Selection pane
--Selection label