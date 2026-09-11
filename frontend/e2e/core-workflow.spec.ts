import AxeBuilder from "@axe-core/playwright";
import { expect, test } from "@playwright/test";

const MULTI_FRAGMENT_INSERT = (
  "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGC"
).repeat(4) + "TAA";

async function createAccount(page, suffix: string) {
  await page.goto("/signup");
  await page.getByLabel("Full name").fill("Dr G Synth");
  await page.getByLabel("Email").fill(`e2e-${suffix}@example.org`);
  await page.getByLabel("Password", { exact: true }).fill("GsynthE2E!2026");
  await page.locator("#password2").fill("GsynthE2E!2026");
  await page.getByRole("button", { name: "Create account" }).click();
  await expect(page.getByRole("heading", { name: /Welcome back/ })).toBeVisible();
}

test("design proceeds through hybridization to restriction cloning", async ({ page }, testInfo) => {
  await createAccount(page, `workflow-${Date.now()}`);
  await page.getByRole("navigation", { name: "Workspace" }).getByRole("link", { name: "Design", exact: true }).click();
  await expect(page.getByRole("heading", { name: "Design a construct" })).toBeVisible();

  await page.getByLabel("Sequence (A/C/G/T)").fill(MULTI_FRAGMENT_INSERT);
  await page.getByRole("button", { name: "Design", exact: true }).click();
  await expect(page.getByText("ESD plan verified", { exact: true })).toBeVisible();
  await expect(page.getByRole("button", { name: /Verify hybridization/ })).toBeDisabled();
  await page.getByRole("button", { name: "Assemble fragments" }).click();
  await expect(page.getByText("Exact reconstruction confirmed", { exact: true })).toBeVisible();
  await expect(page.getByRole("button", { name: "Export all sequences" })).toBeEnabled();
  const sequenceDownload = page.waitForEvent("download");
  await page.getByRole("button", { name: "Export all sequences" }).click();
  expect((await sequenceDownload).suggestedFilename()).toBe("construct_all_sequences.fasta");
  await page.setViewportSize({ width: 320, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await page.getByRole("button", { name: /Verify hybridization/ }).click();

  await expect(page.getByRole("heading", { name: "Hybridization evidence" })).toBeVisible();
  const assembledTop = (await page.locator("#a").inputValue()).replace(/[^A-Za-z]/g, "");
  await page.getByRole("button", { name: /Simulate restriction cloning/ }).click();
  await expect(page).toHaveURL(/\/clone$/);
  await expect(page.getByRole("heading", { name: /Restriction enzyme cloning/ })).toBeVisible();
  await expect(page.getByText(new RegExp(`${assembledTop.length} bp with`))).toBeVisible();
  await expect(page.evaluate(() => window.scrollY)).resolves.toBe(0);
  await page.getByRole("button", { name: "Simulate digestion" }).click();
  await expect(page.getByRole("heading", { name: "Expression reading frame" })).toBeVisible();
  await expect(page.locator(".frame-assessment").getByText("Confirmed", { exact: true })).toBeVisible();
  await expect(page.locator(".frame-assessment").getByText("8 nt to ATG", { exact: true })).toBeVisible();
  await page.getByLabel("Left restriction enzyme", { exact: true }).selectOption("BamHI");
  const configuration = page.getByRole("region", { name: "Cloning configuration" });
  await expect(configuration.getByText("BamHI / XhoI", { exact: true })).toBeVisible();
  await expect(configuration.getByText(/Awaiting simulation/)).toBeVisible();
  await expect(page.getByRole("heading", { name: "Expression reading frame" })).toHaveCount(0);
  await expect(page.getByText(/Cloning at NdeI uses/)).toHaveCount(0);
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await testInfo.attach("updated-cloning-configuration", { body: await configuration.screenshot(), contentType: "image/png" });
});

test("a single-fragment design bypasses assembly and retains its annotated N terminus", async ({ page }, testInfo) => {
  await createAccount(page, `single-fragment-${Date.now()}`);
  await page.getByRole("navigation", { name: "Workspace" })
    .getByRole("link", { name: "Design", exact: true }).click();

  await page.getByRole("button", { name: "Design", exact: true }).click();

  await expect(page.getByText("Single-fragment design verified", { exact: true })).toBeVisible();
  await expect(page.getByRole("heading", { name: "ESD fragment assembly" })).toHaveCount(0);
  await expect(page.getByRole("button", { name: "Assemble fragments" })).toHaveCount(0);
  await expect(page.getByRole("button", { name: "Export all sequences" })).toBeEnabled();
  await expect(page.getByRole("button", { name: /Verify hybridization/ })).toBeEnabled();

  await page.getByRole("button", { name: /Verify hybridization/ }).click();
  await page.getByRole("button", { name: /Simulate restriction cloning/ }).click();
  await page.getByRole("button", { name: "Simulate digestion" }).click();
  await page.getByRole("button", { name: "Ligate compatible ends" }).click();
  await page.getByRole("tab", { name: /^Product/ }).click();
  await page.getByRole("tab", { name: "Sequence", exact: true }).click();

  const annotated = page.getByRole("region", { name: "Coordinate-level annotated sequence" });
  await expect(annotated.getByTitle("Residue 1: Met (ATG)", { exact: true })).toBeVisible();
  for (let residue = 5; residue <= 10; residue += 1) {
    await expect(annotated.getByTitle(`Residue ${residue}: His (CAC)`, { exact: true }).first()).toBeVisible();
  }
  await expect(annotated.getByTitle("Residue 1: Tyr (TAT)", { exact: true })).toHaveCount(0);
  await expect(annotated.getByRole("button", { name: /^construct target, misc_feature/ }).first()).toBeVisible();
  await page.getByRole("button", { name: "Expand viewer", exact: true }).click();
  const expanded = page.getByRole("dialog", { name: "Construct workbench" });
  await expect(expanded).toBeVisible();
  await expanded.getByRole("button", { name: "Whole plasmid" }).click();
  await expanded.getByLabel("Go to base").fill("5400");
  await expanded.getByRole("button", { name: "Go", exact: true }).click();
  await expect(expanded.getByTitle(/^5,400:/)).toBeVisible();
  await expanded.getByRole("button", { name: "Hide details" }).click();
  await expect(expanded.getByRole("complementary", { name: "Selection inspector" })).toHaveCount(0);
  await expect(expanded.getByTitle(/^5,400:/)).toBeInViewport();
  await expanded.getByRole("button", { name: /^Base 5390:/ }).click();
  await expanded.getByRole("button", { name: /^Base 5400:/ }).click({ modifiers: ["Shift"] });
  await expanded.getByRole("button", { name: "Annotate selection" }).click();
  const editor = page.getByRole("dialog", { name: "Annotate a feature" });
  await expect(editor.getByLabel("Start base")).toHaveValue("5390");
  await expect(editor.getByLabel("End base")).toHaveValue("5400");
  await editor.getByLabel("Feature name").fill("User-selected region");
  await editor.getByLabel("Evidence or source note (optional)").fill("Manually selected during the release workflow check.");
  await page.screenshot({ path: testInfo.outputPath("annotation-evidence-editor.png") });
  await expect(editor.getByRole("button", { name: "Add feature" })).toBeInViewport({ ratio: 1 });
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(editor.getByRole("button", { name: "Add feature" })).toBeInViewport({ ratio: 1 });
  await page.screenshot({ path: testInfo.outputPath("annotation-evidence-editor-mobile.png") });
  await page.setViewportSize({ width: 1280, height: 720 });
  await editor.getByRole("button", { name: "Add feature" }).click();
  await expect(expanded.getByRole("button", { name: /^User-selected region, misc_feature/ })).toBeVisible();
  await page.screenshot({ path: testInfo.outputPath("whole-plasmid-expanded.png") });
  await expanded.getByLabel("Go to base").fill("142");
  await expanded.getByRole("button", { name: "Go", exact: true }).click();
  await expanded.getByRole("button", { name: "Show details" }).click();
  const sdMotif = expanded.getByRole("button", { name: /^SD-like motif \(unassigned\), misc_feature, bases 142/ });
  await expect(sdMotif).toHaveText(/\? SD motif/);
  await sdMotif.click();
  await expect(expanded.getByText(/No annotated CDS starts at this codon/)).toBeVisible();
  await expect(expanded.getByText("Detected candidate · function unconfirmed")).toBeVisible();
  await page.screenshot({ path: testInfo.outputPath("sd-context-evidence.png") });
  const results = await new AxeBuilder({ page }).include(".expandable-panel.expanded").analyze();
  expect(results.violations.filter((v) => ["critical", "serious"].includes(v.impact ?? ""))).toEqual([]);
  await page.keyboard.press("Escape");
  await expect(page.getByRole("button", { name: "Expand viewer", exact: true })).toBeFocused();
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);

});

test("mobile navigation is complete and does not widen the viewport", async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await createAccount(page, `mobile-${Date.now()}`);

  const menu = page.getByRole("button", { name: "Menu" });
  await expect(menu).toBeVisible();
  await menu.click();
  const navigation = page.getByRole("navigation", { name: "Workspace" });
  await expect(navigation.getByRole("link", { name: "Restriction clone", exact: true })).toBeVisible();
  await expect(navigation.getByRole("link", { name: "Learn", exact: true })).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
});

test("authentication and home have no serious accessibility violations", async ({ page }) => {
  await page.goto("/login");
  await expect(page.getByText("Designed by Prof. Merzoug Mohamed", { exact: true })).toBeVisible();
  const loginResults = await new AxeBuilder({ page }).analyze();
  expect(loginResults.violations.filter((violation) => ["critical", "serious"].includes(violation.impact ?? ""))).toEqual([]);

  await createAccount(page, `a11y-${Date.now()}`);
  const homeResults = await new AxeBuilder({ page }).analyze();
  expect(homeResults.violations.filter((violation) => ["critical", "serious"].includes(violation.impact ?? ""))).toEqual([]);
});

test("every scientific workspace uses the available width without page overflow", async ({ page }, testInfo) => {
  await createAccount(page, `layouts-${Date.now()}`);
  const workspaces = [
    "/optimise", "/design", "/pcr", "/hybridize", "/clone", "/verify", "/align",
  ];

  for (const route of workspaces) {
    await page.goto(route);
    const layout = page.locator(".design-layout").first();
    await expect(layout).toBeVisible();
    expect(await layout.evaluate((element) => getComputedStyle(element).gridTemplateColumns.split(" ").length)).toBe(1);
    await page.screenshot({ path: testInfo.outputPath(`${route.slice(1)}-desktop.png`), fullPage: true });
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  }

  for (const width of [768, 390, 320]) {
    await page.setViewportSize({ width, height: 844 });
    for (const route of workspaces) {
      await page.goto(route);
      await expect(page.locator(".design-layout").first()).toBeVisible();
      await page.screenshot({ path: testInfo.outputPath(`${route.slice(1)}-${width}.png`), fullPage: true });
      expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
    }
  }

  await page.goto("/hybridize");
  const stageLabels = page.locator(".core-workflow-stage strong");
  await expect(stageLabels).toHaveCount(3);
  expect(await stageLabels.evaluateAll((labels) => labels.every(
    (label) => label.scrollWidth <= label.clientWidth && label.scrollHeight <= label.clientHeight,
  ))).toBe(true);
});

test("edited cloning primers are revalidated before cloning", async ({ page }) => {
  await createAccount(page, `primers-${Date.now()}`);
  await page.getByRole("navigation", { name: "Workspace" })
    .getByRole("link", { name: "Primers & PCR", exact: true }).click();
  await page.getByLabel("Sequence to amplify (A/C/G/T)").fill(
    "ATGACCACCAGCAAACTGGGCAAAGGCCTGGGCTATATTGGCAACAACGGCGCGCACATGGGCTTAAACTTAGCATTACTGGGCCTGGCGAGCCTGCTGGGCAAAGGCATTAGCAAACTGGGC",
  );
  await page.getByRole("button", { name: "Design PCR" }).click();
  await expect(page.getByText("G-Synth design", { exact: true })).toBeVisible();

  await page.getByRole("button", { name: "Edit primers" }).click();
  await page.getByLabel("Forward primer (5′→3′)").fill(
    "AACCGGCATATGACCACCAGCAAACTGGGC",
  );
  await expect(page.getByText("Edited primers need validation")).toBeVisible();
  await page.getByRole("button", { name: "Revalidate primers" }).click();

  await expect(page.getByText("Edited · validated", { exact: true })).toBeVisible();
  await expect(page.getByLabel("Nucleotide-level digested insert and cohesive ends")).toBeVisible();
  await expect(page.getByText("5′-TA", { exact: true })).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);

  const results = await new AxeBuilder({ page }).analyze();
  expect(results.violations.filter((violation) => ["critical", "serious"].includes(violation.impact ?? ""))).toEqual([]);
});

test("workspace text remains readable when enlarged", async ({ page }, testInfo) => {
  await createAccount(page, `readability-${Date.now()}`);
  await page.setViewportSize({ width: 1280, height: 900 });
  for (const route of ["/", "/learn", "/design", "/clone", "/verify", "/align"]) {
    await page.goto(route);
    await expect(page.locator("main")).toBeVisible();
    expect(await page.locator(".note, .field-hint, .field > label, .card small").evaluateAll((elements) =>
      elements.filter((element) => element.getClientRects().length > 0).every((element) =>
        parseFloat(getComputedStyle(element).fontSize) >= 14,
      ),
    )).toBe(true);
    const enlarged = await page.addStyleTag({ content: "html { font-size: 200%; }" });
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
    expect(await page.locator(".pipeline, .rail-group-label, .core-workflow-trail").evaluateAll((elements) =>
      elements.filter((element) => element.getClientRects().length > 0).every((element) => element.scrollWidth <= element.clientWidth),
    )).toBe(true);
    await page.screenshot({ path: testInfo.outputPath(`${route.slice(1) || "home"}-text-200.png`), fullPage: true });
    await enlarged.evaluate((element) => element.remove());
  }
});
