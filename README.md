> **Research prototype:** [`enzyme-evidence/`](enzyme-evidence/) — evidence integration for enzyme EC classification on the CARE benchmark (retrieval + learned integrators, leakage-safe evaluation).

Projektidee: BioCloud – Eine KI-gestützte Bioinformatik-Plattform
Ziel
Entwicklung einer Web-App mit Flask, die FASTA-Dateien hochlädt, Sequenzen analysiert, visualisiert, klassifiziert und mit synthetischen Daten trainiert, während die gesamte Pipeline auf AWS mit Docker läuft.

Schrittweise Umsetzung
Hier sind die wichtigsten Schritte, die dir helfen, verschiedene Technologien zu erlernen.

1️⃣ Flask Web App (Frontend & API)
Flask-App mit Datei-Upload für FASTA- und FASTQ-Sequenzen.
REST API, die Daten an AWS sendet.
Visualisierung der Sequenzdaten als Multiple Sequence Alignment (MSA), Phylogenetischer Baum oder 3D-Proteinstruktur mit Matplotlib, Biopython oder Plotly.
✅ Tech: Flask, HTML/CSS, Bootstrap, Matplotlib, Plotly

2️⃣ AWS Backend (Python-Skript & Bioinformatik-Pipeline)
AWS Lambda oder AWS Batch mit EC2, um FASTA-Dateien in die Cloud zu laden.
Alignment mit Clustal Omega oder MAFFT in einem Bash-Skript auf AWS.
Python-Skript für weitere Bioinformatik-Analysen, z. B. Proteinstruktur-Vorhersage mit AlphaFold oder Sekundärstruktur-Vorhersage mit biopython.
✅ Tech: AWS Lambda, S3, Bash, Biopython, AWS Batch

3️⃣ Machine Learning für Klassifikation (ML-Modell in Python)
Trainiere ein ML-Modell zur Klassifikation von DNA/Protein-Sequenzen.
Nutze scikit-learn oder TensorFlow für die Klassifikation.
Modell trainieren auf realen und synthetischen Daten (siehe Schritt 5).
✅ Tech: Pandas, Scikit-Learn, TensorFlow, Numpy

4️⃣ Deep Learning für Sequenz-Vorhersage (AI-Modelle nutzen)
Nutze einen Transformer oder CNN, um DNA-Sequenzen zu analysieren.
Vergleich von CNNs, RNNs und Transformers für DNA-Klassifikation.
Trainiere dein Modell auf AWS mit GPU-Unterstützung (z. B. SageMaker oder EC2 mit CUDA).
✅ Tech: TensorFlow/Keras, PyTorch, AWS SageMaker

5️⃣ AI-Synthese von DNA/Protein-Sequenzen für Training (AI Synth Data)
Erzeuge synthetische DNA-Sequenzen mit Generative Adversarial Networks (GANs) oder LSTMs.
Vergleiche generierte Sequenzen mit realen Daten und evaluiere die Ähnlichkeit.
✅ Tech: PyTorch, TensorFlow, Generative AI (GANs, VAEs)

6️⃣ RAG (Retrieval-Augmented Generation) für Bioinformatik-Daten
Bau eine RAG (z. B. mit LangChain), die auf Bioinformatik-Daten spezialisiert ist.
Nutze FAISS für Vektorsuche in DNA/Protein-Datenbanken.
Ermögliche eine interaktive Abfrage biologischer Daten mit LLMs.
✅ Tech: LangChain, FAISS, OpenAI API, HuggingFace Transformers

7️⃣ Docker & Deployment auf AWS
Erstelle ein Docker-Image der gesamten App (Flask + ML + AWS Integration).
Nutze AWS Fargate oder Kubernetes für Skalierung.
Automatisiere das Deployment mit CI/CD (GitHub Actions oder AWS CodePipeline).
✅ Tech: Docker, Kubernetes, AWS ECS, CI/CD (GitHub Actions, AWS CodeBuild)

Endergebnis
✨ Eine Cloud-gestützte Bioinformatik-Plattform, die:

Sequenzen hochlädt und analysiert (AWS Lambda + Bash).
ML-Modelle zur Klassifikation nutzt.
Deep Learning für Sequenzvorhersagen einsetzt.
AI-generierte Sequenzen erstellt (SynthData).
RAG nutzt, um Bioinformatik-Daten zu durchsuchen.
Als Docker-Container läuft und auf AWS deployt wird.
