#!/bin/bash
# Launcher for Diagram Viewer HTML page

echo "==================================================="
echo "  Diagram Viewer Launcher"
echo "==================================================="
echo ""

# Check if we're in WSL
if grep -qi microsoft /proc/version; then
    echo "Detected WSL environment"
    
    # Get the absolute path and convert to Windows format
    HTML_PATH="$(realpath /home/matthematics/schubmult/diagram_viewer.html)"
    WINDOWS_PATH="$(wslpath -w "$HTML_PATH")"
    
    echo "Opening in Windows browser..."
    # Use powershell.exe to open the file (works better than cmd.exe with UNC paths)
    powershell.exe -Command "Start-Process '$WINDOWS_PATH'"
    
    echo ""
    echo "Alternative methods:"
    echo "-------------------"
    echo "1. Run Python GUI version:"
    echo "   python little_bump_viewer.py"
    echo ""
    echo "2. Run original Java applet (if you have JDK):"
    echo "   cd /mnt/c/users/matth/downloads/diagramviewer"
    echo "   appletviewer diagram_viewer.html"
    echo ""
    echo "3. Run as standalone JAR:"
    echo "   cd /mnt/c/users/matth/downloads/diagramviewer"
    echo "   java -jar DiagramViewer.jar"
    
else
    echo "Opening in default browser..."
    xdg-open "/home/matthematics/schubmult/diagram_viewer.html"
fi

echo ""
echo "==================================================="
