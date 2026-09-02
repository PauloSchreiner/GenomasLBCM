import sys
import os
import pandas as pd

def main():
    tsv_file = sys.argv[1]
    out_html = sys.argv[2]
    txt_files = sys.argv[3:]

    # Loads run_history tsv
    try:
        df = pd.read_csv(tsv_file, sep='\t', keep_default_na=False)
        html_table = df.to_html(index=False, classes='overview-table', border=0)
    except Exception as e:
        html_table = f"<p>Error loading .tsv file {e}</p>"

    # Prepares tabs
    tab_buttons = '<button class="tablinks active" onclick="openTab(event, \'General\')">Overview</button>\n'
    tab_contents = f'<div id="General" class="tabcontent" style="display:block;"><h2>Run history</h2><div style="overflow-x:auto;">{html_table}</div></div>\n'

    for txt in txt_files:
        sample_name = os.path.basename(txt).replace('_assembly_summary.txt', '')
        try:
            with open(txt, 'r', encoding='utf-8') as f:
                content = f.read()
        except Exception as e:
            content = f"Erro ao ler arquivo: {e}"
            
        tab_buttons += f'<button class="tablinks" onclick="openTab(event, \'{sample_name}\')">{sample_name}</button>\n'
        tab_contents += f'<div id="{sample_name}" class="tabcontent"><h2>{sample_name}</h2><pre class="txt-content">{content}</pre></div>\n'

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Assembly report</title>
<style>
    body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #f4f6f9; color: #333; margin: 0; padding: 20px; }}
    h1 {{ color: #2c3e50; text-align: center; }}
    .tab {{ overflow: hidden; border: 1px solid #ccc; background-color: #e9ecef; border-radius: 8px 8px 0 0; display: flex; flex-wrap: wrap; }}
    .tab button {{ background-color: inherit; border: none; outline: none; cursor: pointer; padding: 12px 18px; transition: 0.3s; font-size: 14px; font-weight: bold; color: #495057; }}
    .tab button:hover {{ background-color: #dee2e6; }}
    .tab button.active {{ background-color: #fff; color: #007bff; border-bottom: 3px solid #007bff; }}
    .tabcontent {{ display: none; padding: 20px; border: 1px solid #ccc; border-top: none; background-color: #fff; border-radius: 0 0 8px 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); min-height: 500px; }}
    .overview-table {{ width: 100%; border-collapse: collapse; margin-top: 15px; font-size: 13px; white-space: nowrap; }}
    .overview-table th {{ background-color: #007bff; color: white; padding: 10px; text-align: left; position: sticky; top: 0; }}
    .overview-table td {{ border-bottom: 1px solid #ddd; padding: 8px; }}
    .overview-table tr:nth-child(even) {{ background-color: #f8f9fa; }}
    .overview-table tr:hover {{ background-color: #e2e6ea; }}
    .txt-content {{ background-color: #1e1e1e; color: #d4d4d4; padding: 18px; border-radius: 6px; overflow-x: auto; font-family: 'Courier New', Courier, monospace; font-size: 13px; line-height: 1.4; }}
</style>
</head>
<body>
<h1>GenomasLBCM - Assembly report</h1>
<div class="tab">
  {tab_buttons}
</div>
{tab_contents}
<script>
function openTab(evt, tabName) {{
  var i, tabcontent, tablinks;
  tabcontent = document.getElementsByClassName("tabcontent");
  for (i = 0; i < tabcontent.length; i++) {{ tabcontent[i].style.display = "none"; }}
  tablinks = document.getElementsByClassName("tablinks");
  for (i = 0; i < tablinks.length; i++) {{ tablinks[i].className = tablinks[i].className.replace(" active", ""); }}
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}}
</script>
</body>
</html>
"""
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html_template)

if __name__ == "__main__":
    main()