import sys
import os
import sqlite3
import snapgene_reader
from Bio import SeqIO
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QLabel, QPushButton, QVBoxLayout, QWidget,
    QFileDialog, QMessageBox, QDialog, QLineEdit, QTextEdit, QInputDialog
)
from PyQt6.QtCore import Qt

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPixmap

# dna_features_viewer for circular/plasmid map visualization
from dna_features_viewer import GraphicRecord, GraphicFeature

DB_PATH = "projects.db"  # Database file


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("DNA Cloning Tracker")
        self.setGeometry(100, 100, 800, 600)

        self.current_project = None
        self.current_sequence = None

        # Create the main layout
        self.layout = QVBoxLayout()

        # Title label
        self.label = QLabel("Welcome to DNA Cloning Tracker", self)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(self.label)

        # Create New Project button
        self.new_project_btn = QPushButton("Create New Project", self)
        self.new_project_btn.clicked.connect(self.create_new_project)
        self.layout.addWidget(self.new_project_btn)

        # Upload DNA Map button
        self.upload_btn = QPushButton("Upload DNA Map", self)
        self.upload_btn.clicked.connect(self.open_dna_file)
        self.upload_btn.setEnabled(False)
        self.layout.addWidget(self.upload_btn)

        # Text box to display basic sequence info
        self.dna_details = QTextEdit(self)
        self.dna_details.setReadOnly(True)
        self.layout.addWidget(self.dna_details)

        # **QLabel for displaying the circular plasmid image**
        self.plasmid_label = QLabel(self)
        self.layout.addWidget(self.plasmid_label)

        # Optionally, other buttons (View Gel, Track Sequences, Notes)...
        # self.view_gel_btn = ...
        # self.track_sequences_btn = ...
        # self.notes_btn = ...

        # Finalize layout
        central_widget = QWidget(self)
        central_widget.setLayout(self.layout)
        self.setCentralWidget(central_widget)

    def init_database(self):
        """Initialize SQLite database for projects."""
        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS projects (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            notes TEXT
        )""")
        conn.commit()
        conn.close()

    def create_new_project(self):
        """Creates a new project by setting its name."""
        project_name, ok = QInputDialog.getText(self, "Create New Project",
                                                "Enter project name:")
        if ok and project_name.strip():
            self.current_project = project_name.strip()
            self.label.setText(f"Project: {self.current_project}")
            self.upload_btn.setEnabled(True)  # Enable file upload
            QMessageBox.information(self, "Success",
                                    f"Project '{self.current_project}' created!")
        else:
            QMessageBox.warning(self, "Error", "Project name cannot be empty.")

    def save_project(self, dialog, project_name):
        """Saves the new project in the database."""
        if not project_name.strip():
            QMessageBox.warning(self, "Error", "Project name cannot be empty.")
            return

        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        try:
            cursor.execute("INSERT INTO projects (name, notes) VALUES (?, ?)",
                           (project_name, ""))
            conn.commit()
            self.current_project = project_name
            self.label.setText(f"Project: {project_name}")
            self.enable_project_features()
            QMessageBox.information(self, "Success",
                                    f"Project '{project_name}' created!")
            dialog.accept()
        except sqlite3.IntegrityError:
            QMessageBox.warning(self, "Error",
                                "A project with this name already exists.")
        finally:
            conn.close()

    def enable_project_features(self):
        """Enable features once a project is created."""
        self.upload_btn.setEnabled(True)
        self.view_gel_btn.setEnabled(True)
        self.track_sequences_btn.setEnabled(True)
        self.notes_btn.setEnabled(True)

    def open_dna_file(self):
        """Opens and processes a DNA sequence file after creating a project."""
        if not self.current_project:
            QMessageBox.warning(self, "Error",
                                "Please create a project first.")
            return

        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self,
            "Open DNA File",
            "",
            "SnapGene Files (*.dna);;GenBank Files (*.gb);;FASTA Files (*.fasta);;All Files (*)"
        )

        if file_path:
            self.label.setText(
                f"Project: {self.current_project}\nLoaded: {os.path.basename(file_path)}")
            self.process_dna_file(file_path)
        else:
            self.label.setText(
                f"Project: {self.current_project}\nNo file selected.")

    def process_dna_file(self, file_path):
        """Reads and extracts DNA sequence information from .dna, .gb, and .fasta files."""
        try:
            if file_path.endswith(".gb"):
                record = SeqIO.read(file_path, "genbank")
                # Handle GenBank as before...
                ...
            elif file_path.endswith(".fasta"):
                record = SeqIO.read(file_path, "fasta")
                # Handle FASTA as before...
                ...
            elif file_path.endswith(".dna"):
                data = snapgene_reader.snapgene_file_to_dict(file_path)
                sequence = data["seq"]
                name = data.get("name", "Unknown")

                # Display basic info in the text box
                sequence_text = (
                    f"Name: {name}\n"
                    f"Length: {len(sequence)} bp\n\n"
                    f"Sequence:\n{sequence[:100]}..."
                )
                self.dna_details.setText(sequence_text)

                # Draw the circular plasmid map
                self.show_circular_map(data)
                return
            else:
                QMessageBox.warning(
                    self, "Unsupported Format",
                    "Only .dna, .gb, and .fasta files are supported."
                )
                return

            # If it's .gb or .fasta, handle as you already do...
            ...

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"Failed to process file:\n{e}")

    def show_circular_map(self, data):
        """Use dna_features_viewer to draw a circular plasmid map from SnapGene data."""
        # 1) Build a list of GraphicFeature objects from SnapGene 'features'
        features = []
        for feat in data.get("features", []):
            start = feat["start"]
            end = feat["end"]
            # SnapGene might store direction as "forward" or "reverse"
            direction = feat.get("direction", "forward").lower()
            strand = +1 if direction == "forward" else -1

            # Use the 'type' or 'name' as the label
            label = feat.get("type", "Feature")

            # Create a dna_features_viewer GraphicFeature
            features.append(
                GraphicFeature(start=start, end=end, strand=strand,
                               label=label)
            )

        # 2) Create a GraphicRecord for the circular plasmid
        record = GraphicRecord(
            sequence_length=len(data["seq"]),
            features=features
        )

        # 3) Plot a circular figure
        circular_diagram = record.plot_circular(
            figure_width=6)  # Adjust figure_width as needed
        circular_diagram.figure.tight_layout()

        # 4) Save the figure to a temporary file (overwrite each time)
        output_path = "plasmid.png"
        circular_diagram.figure.savefig(output_path, bbox_inches="tight")
        circular_diagram.figure.clf()  # Clear the figure from memory

        # 5) Load that image into the QLabel
        pixmap = QPixmap(output_path)
        self.plasmid_label.setPixmap(pixmap)

    def open_notes(self):
        """Opens a Notion-style note editor."""
        if not self.current_project:
            QMessageBox.warning(self, "Error", "No project selected.")
            return

        notes_dialog = NotesWindow(self.current_project)
        notes_dialog.exec()

    def show_message(self):
        """Placeholder for future features."""
        QMessageBox.information(self, "Coming Soon",
                                "This feature is under development!")


class NotesWindow(QDialog):
    """A window for taking Notion-style project notes."""

    def __init__(self, project_name):
        super().__init__()

        self.project_name = project_name
        self.setWindowTitle(f"Notes: {project_name}")
        self.setGeometry(300, 200, 600, 400)

        layout = QVBoxLayout()

        self.notes_field = QTextEdit(self)
        self.notes_field.setPlaceholderText(
            "Write your project notes, protocols, or observations here...")
        layout.addWidget(self.notes_field)

        save_btn = QPushButton("Save Notes", self)
        save_btn.clicked.connect(self.save_notes)
        layout.addWidget(save_btn)

        self.setLayout(layout)
        self.load_existing_notes()

    def load_existing_notes(self):
        """Loads notes from the database."""
        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute("SELECT notes FROM projects WHERE name = ?",
                       (self.project_name,))
        notes = cursor.fetchone()
        conn.close()

        if notes and notes[0]:
            self.notes_field.setText(notes[0])

    def save_notes(self):
        """Saves notes to the database."""
        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute("UPDATE projects SET notes = ? WHERE name = ?",
                       (self.notes_field.toPlainText(), self.project_name))
        conn.commit()
        conn.close()
        QMessageBox.information(self, "Saved", "Notes saved successfully!")


# Run the app
app = QApplication(sys.argv)
window = MainWindow()
window.show()
sys.exit(app.exec())
