"""Append four slides on the ovary foundation-model project to the lab-update deck, in its style."""
import sys
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

src, dst = sys.argv[1], sys.argv[2]
prs = Presentation(src)
LAYOUT = prs.slide_layouts[1]           # 'Title and body'
BG, CARD, CARD_LINE = "F7F6F2", "FDFDFB", "E3E1DA"
INK, MUTED, FOOT, GREY = "1A1A1A", "4A5568", "6A7179", "898781"
BLUE, ORANGE, LIGHT_BLUE, SLATE = "2A78D6", "EB6834", "9EC5F4", "515150"
FOOTER = "Visium-QC · lab meeting, November 2026 · {n}"


def rgb(h):
    return RGBColor.from_string(h)


def new_slide(title, n):
    s = prs.slides.add_slide(LAYOUT)
    s.background.fill.solid(); s.background.fill.fore_color.rgb = rgb(BG)
    for ph in list(s.placeholders):
        if ph.placeholder_format.idx != 0:
            ph._element.getparent().remove(ph._element)
    t = s.shapes.title
    t.left, t.top, t.width, t.height = Inches(0.89), Inches(0.89), Inches(11.9), Inches(0.96)
    tf = t.text_frame; tf.text = title
    for p in tf.paragraphs:
        for r in p.runs:
            r.font.name = "Georgia"; r.font.size = Pt(30); r.font.bold = True; r.font.color.rgb = rgb(INK)
    text(s, 0.89, 6.82, 4.4, 0.26, FOOTER.format(n=n), size=12, color=FOOT)
    return s


def text(s, x, y, w, h, content, size=13, color=INK, bold=False, font="Arial", align=PP_ALIGN.LEFT,
         anchor=MSO_ANCHOR.TOP, space_after=6, name=None):
    """content: a string, or a list of paragraphs; a paragraph is a string or a list of (text, {opts}) runs."""
    box = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        box.name = name
    tf = box.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.05); tf.margin_top = tf.margin_bottom = Inches(0.03)
    tf.vertical_anchor = anchor
    paras = content if isinstance(content, list) else [content]
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align; p.space_after = Pt(space_after)
        runs = para if isinstance(para, list) else [(para, {})]
        for txt, opts in runs:
            r = p.add_run(); r.text = txt
            r.font.name = opts.get("font", font); r.font.size = Pt(opts.get("size", size))
            r.font.bold = opts.get("bold", bold); r.font.italic = opts.get("italic", False)
            r.font.color.rgb = rgb(opts.get("color", color))
    return box


def card(s, x, y, w, h, fill=CARD, line=CARD_LINE):
    c = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    c.adjustments[0] = 0.035
    c.fill.solid(); c.fill.fore_color.rgb = rgb(fill)
    c.line.color.rgb = rgb(line); c.line.width = Pt(0.75)
    c.shadow.inherit = False
    c.text_frame.text = ""
    return c


def arrow(s, x, y, w, h, fill=LIGHT_BLUE):
    a = s.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(x), Inches(y), Inches(w), Inches(h))
    a.fill.solid(); a.fill.fore_color.rgb = rgb(fill); a.line.fill.background(); a.shadow.inherit = False
    return a


def box(s, x, y, w, h, label, body, fill, label_color="FFFFFF", body_color="FFFFFF", size=12.5):
    c = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    c.adjustments[0] = 0.06
    c.fill.solid(); c.fill.fore_color.rgb = rgb(fill); c.line.fill.background(); c.shadow.inherit = False
    tf = c.text_frame; tf.word_wrap = True; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = Inches(0.12)
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = label; r.font.name = "Arial"; r.font.size = Pt(12); r.font.bold = True; r.font.color.rgb = rgb(label_color)
    if body:
        p2 = tf.add_paragraph(); p2.alignment = PP_ALIGN.CENTER
        r2 = p2.add_run(); r2.text = body; r2.font.name = "Arial"; r2.font.size = Pt(size); r2.font.color.rgb = rgb(body_color)
    return c


n0 = len(prs.slides)

# ---------------------------------------------------------------- slide A: why
s = new_slide("Next: a model that knows the ovary's cells", n0 + 1)
stats = [("0.4%", "of the public cell atlas that today's models learn from is human ovary, and half of that is tumour"),
         ("0.08%", "is mouse ovary: two studies, no nuclei, no estrous cycle, no ovary older than 9 months"),
         ("237", "mouse oocytes in the whole atlas. The follicle, the part that defines the organ, is almost absent")]
for i, (num, cap) in enumerate(stats):
    x = 0.89 + i * 3.93
    card(s, x, 2.05, 3.7, 2.15)
    text(s, x + 0.28, 2.3, 3.23, 0.86, num, size=60, color=BLUE, bold=True)
    text(s, x + 0.28, 3.27, 3.23, 0.9, cap, size=13, color=MUTED)
text(s, 0.89, 4.5, 5.65, 1.9, [
    [("What a foundation model is. ", {"bold": True}),
     ("A model trained on millions of single cells from public atlases. It learns what each cell type and state looks like, "
      "so it can label new data and, in the best case, predict how cells change under a treatment. The Arc Institute now runs a "
      "yearly contest for exactly that prediction.", {})]], size=13.5, color=INK)
text(s, 6.81, 4.5, 5.63, 1.9, [
    [("Why the ovary needs its own. ", {"bold": True}),
     ("Those models have seen almost no ovary, and none of the biology we ask about: the cycle, ovulation, ageing, the response to "
      "surgery. The data exist in GEO, more than a million mouse cells, but nobody has gathered them into something a model can learn from.", {})]],
     size=13.5, color=INK)
s.notes_slide.notes_text_frame.text = (
    "Framing for a non-computational audience: a foundation model is a model pre-trained on millions of cells that can be reused for "
    "many tasks. The numbers come from the CELLxGENE Census release of November 2025, queried directly. The Arc Institute's Virtual "
    "Cell Challenge (2026 edition) is the motivation and the yardstick; we are not entering it.")

# ---------------------------------------------------------------- slide B: what
s = new_slide("What we will build: three things, each useful alone", n0 + 2)
cols = [
    ("1 · THE COLLECTION", BLUE,
     "Every public mouse ovary dataset we can find, re-processed from raw counts into one collection with one shared vocabulary of "
     "cell types: the estrous cycle, ovulation time courses after hCG, ageing and menopause models, and our own oophorectomy nuclei. "
     "About one to two million cells. Most of it sits in GEO and never reached the atlases the big models learn from."),
    ("2 · THE TEST", ORANGE,
     "A held-out test. The model is given only a compartment's baseline cells, say granulosa before hCG, or young stroma, and asked "
     "to predict the treated state: 4 hours after hCG, or aged. Its answer is scored against the real cells on a 0-to-1 scale, where "
     "0 is no change at all and 1 is what a second animal achieves."),
    ("3 · THE MODEL", BLUE,
     "A modest mouse model, started from a public one and trained further on the ovary collection, judged against simple methods "
     "on the same test. If it does not beat them by a margin we have written down in advance, the collection and the test are "
     "still the papers, and we say so."),
]
for i, (label, col, body) in enumerate(cols):
    x = 0.89 + i * 3.93
    card(s, x, 2.05, 3.7, 3.1)
    text(s, x + 0.28, 2.3, 3.23, 0.26, label, size=12, color=col, bold=True)
    text(s, x + 0.28, 2.68, 3.23, 2.4, body, size=13, color=INK)
text(s, 0.89, 5.5, 11.44, 0.7,
     "The yardstick is the Arc Institute's 2026 Virtual Cell Challenge: the same task shape and the same way of scoring, "
     "copied onto ovarian compartments and conditions. We use its rules; we are not entering it.", size=14, color=MUTED)
s.notes_slide.notes_text_frame.text = (
    "Datasets already in hand: Isola 2024 (ageing, Stout lab), Huang 2025 (ovulation, PLOS Biology), Wang and Suh 2026 (estropause, PNAS), "
    "the Schuh lab ovulation series (Nature Cell Biology 2024), Gaylord 2025 (Science), Morris 2022 (the cycling ovary, Pépin lab, pending "
    "a download), plus our 30-day oophorectomy nuclei. The base model is scPRINT-2 (Apache-2.0 weights); the primary comparison is scVI/scANVI.")

# ---------------------------------------------------------------- slide C: how the test works
s = new_slide("How the test works, in one picture", n0 + 3)
y0 = 2.15; h = 1.25
box(s, 0.89, y0, 2.55, h, "BASELINE CELLS", "one compartment, before the treatment (granulosa at 0 h; young stroma)", "DDE9F8", BLUE, INK)
arrow(s, 3.52, y0 + 0.42, 0.55, 0.42)
box(s, 4.15, y0, 2.2, h, "THE MODEL", "sees only those cells and the name of the treatment", BLUE)
arrow(s, 6.43, y0 + 0.42, 0.55, 0.42)
box(s, 7.06, y0, 2.55, h, "PREDICTED CELLS", "what it thinks 4 h after hCG, or aged, looks like", "DDE9F8", BLUE, INK)
text(s, 9.66, y0 + 0.3, 0.6, 0.6, "vs", size=18, color=MUTED, bold=True, align=PP_ALIGN.CENTER)
box(s, 10.27, y0, 2.2, h, "REAL CELLS", "the measured treated cells, held back from the model", "FBE3D9", ORANGE, INK)
# the scale
text(s, 0.89, 3.75, 11.6, 0.3, "The score, per compartment and treatment", size=13, color=GREY, bold=True)
line = s.shapes.add_connector(1, Inches(1.4), Inches(4.45), Inches(11.95), Inches(4.45))
line.line.color.rgb = rgb(SLATE); line.line.width = Pt(2.25)
for xpos, colour, big, small, align, cx, cw in [
    (1.4, SLATE, "0", "the baseline copied: no change predicted", PP_ALIGN.LEFT, 0.89, 3.4),
    (7.0, ORANGE, "the bar", "what you get by borrowing the response from another dataset, no learning at all", PP_ALIGN.CENTER, 5.3, 3.4),
    (11.95, BLUE, "1", "what a second animal reproduces; the ceiling we can measure", PP_ALIGN.RIGHT, 9.04, 3.4),
]:
    dot = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(xpos - 0.12), Inches(4.33), Inches(0.24), Inches(0.24))
    dot.fill.solid(); dot.fill.fore_color.rgb = rgb(colour); dot.line.fill.background(); dot.shadow.inherit = False
    text(s, xpos - 0.6, 3.98, 1.2, 0.32, big, size=16, color=colour, bold=True, align=PP_ALIGN.CENTER)
    text(s, cx, 4.68, cw, 0.75, small, size=12, color=MUTED, align=align)
text(s, 0.89, 5.5, 11.44, 1.15, [
    "Whole datasets are held out, never one condition of a dataset. The genes that name the treatment (the immediate-early genes after hCG) are "
    "removed from the scoring, as the Challenge removes the knocked-down gene. A compartment that does not respond is listed but not scored.",
    [("The rule, written down before any model is trained: ", {"bold": True}),
     ("a model counts only if it clears the bar in at least two held-out contexts. Otherwise we publish the collection, the test and the ceilings, "
      "and report the model as the negative result it is.", {})]], size=13, color=INK)
s.notes_slide.notes_text_frame.text = (
    "The scoring copies the Challenge's 2026 contract: six measures of how well the predicted cells match the real ones (which genes change, "
    "in which direction, by how much; how separable the responses are), each rescaled so that the control copy is 0 and a split-half replicate "
    "across animals is 1. The bar is retrieval and transfer, the approach that placed in both Challenges. Preregistration lives in the repo.")

# ---------------------------------------------------------------- slide D: our data, status, asks
s = new_slide("Our oophorectomy data, status, and the asks", n0 + 4)
cols = [
    ("OUR DATA", BLUE,
     "The 30-day ooph vs sham nuclei become a test the model never trains on: no public single-cell dataset of compensatory ovulation "
     "exists, so ours is the only one. The 120-day cohort is designed to be scorable: at least 400 nuclei per compartment per arm, three "
     "mice per arm if the budget allows, estrous stage recorded at collection."),
    ("NOW", ORANGE,
     "The pipeline is built and tested. Four public datasets are through QC and labelling, and the first ceilings table ran: 45 "
     "compartment-by-condition tests. It taught us two things: one deposit is raw droplets and now gets cell calling; and with "
     "thousands of cells per group every gene looks significant, so groups are now capped, as the Challenge does."),
    ("ASKS", BLUE,
     "One hour on the cell-type crosswalk: our names against the Cell Ontology's new granulosa, theca and luteal terms, with the markers "
     "beside each. A decision on the 120-day design. And a conversation about the one thing no lab has: a perturbation experiment in the "
     "ovary itself, which is what moved the score in every Challenge so far."),
]
for i, (label, col, body) in enumerate(cols):
    x = 0.89 + i * 3.93
    card(s, x, 2.05, 3.7, 3.1)
    text(s, x + 0.28, 2.3, 3.23, 0.26, label, size=12, color=col, bold=True)
    text(s, x + 0.28, 2.68, 3.23, 2.4, body, size=13, color=INK)
text(s, 0.89, 5.5, 11.44, 0.7,
     "Timeline: the last two atlases and the corrected tables this month; the collection and the first results, November to January; "
     "the model and its ablations, spring; the write-up, next autumn. Everything is on Biowulf and in the Virtual_Cell repository.",
     size=14, color=MUTED)
s.notes_slide.notes_text_frame.text = (
    "Status on 2026-11-05: Huang 2025 (25,080 cells, 16 named clusters), Wang and Suh 2026 (78,538 cells, 18 ovaries), Gaylord 2025 "
    "(51,841 cells, 6 mice) and the Schuh 2024 follicle series are through QC; Isola 2024 (Census) and Morris 2022 (OSF) are pending. "
    "The first ceilings run (45 bundles, two hours) showed the Schuh deposit holds every droplet (49 million barcodes) and that uncapped "
    "groups make every gene significant; both are fixed in the pipeline and the table is being rerun. "
    "The crosswalk review is the one task that needs the lab's biology rather than compute. The in vivo perturbation point is an option, not a "
    "dependency: the plan stands without it, but every 2026 result says in-context perturbation data is what makes a learned model beat transfer.")

prs.save(dst)
print(f"wrote {dst}: {len(prs.slides)} slides")
