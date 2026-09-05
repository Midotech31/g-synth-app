import AxeBuilder from "@axe-core/playwright";
import { expect, test } from "@playwright/test";

async function createAccount(page, suffix: string) {
  await page.goto("/signup");
  await page.getByLabel("Full name").fill("Dr G Synth");
  await page.getByLabel("Email").fill(`e2e-${suffix}@example.org`);
  await page.getByLabel("Password", { exact: true }).fill("GsynthE2E!2026");
  await page.locator("#password2").fill("GsynthE2E!2026");
  await page.getByRole("button", { name: "Create account" }).click();
  await expect(page.getByRole("heading", { name: /Welcome back/ })).toBeVisible();
}

test("design proceeds through hybridization to restriction cloning", async ({ page }) => {
  await createAccount(page, `workflow-${Date.now()}`);
  await page.getByRole("navigation", { name: "Workspace" }).getByRole("link", { name: "Design", exact: true }).click();
  await expect(page.getByRole("heading", { name: "Design a construct" })).toBeVisible();

  await page.getByRole("button", { name: "Design", exact: true }).click();
  await expect(page.getByText("Assembly verified", { exact: true })).toBeVisible();
  await page.getByRole("button", { name: /Verify hybridization/ }).click();

  await expect(page.getByRole("heading", { name: "Hybridization evidence" })).toBeVisible();
  await page.getByRole("button", { name: /Simulate restriction cloning/ }).click();
  await expect(page).toHaveURL(/\/clone$/);
  await expect(page.getByRole("heading", { name: /Restriction enzyme cloning/ })).toBeVisible();
  await expect(page.evaluate(() => window.scrollY)).resolves.toBe(0);
  await page.getByRole("button", { name: "Simulate digestion" }).click();
  await expect(page.getByRole("heading", { name: "Expression reading frame" })).toBeVisible();
  await expect(page.locator(".frame-assessment").getByText("Confirmed", { exact: true })).toBeVisible();
  await expect(page.locator(".frame-assessment").getByText("8 nt to ATG", { exact: true })).toBeVisible();
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

test("every scientific workspace uses the available width without page overflow", async ({ page }) => {
  await createAccount(page, `layouts-${Date.now()}`);
  const workspaces = [
    "/optimise", "/design", "/pcr", "/hybridize", "/clone", "/verify", "/align",
  ];

  for (const route of workspaces) {
    await page.goto(route);
    const layout = page.locator(".design-layout").first();
    await expect(layout).toBeVisible();
    expect(await layout.evaluate((element) => getComputedStyle(element).gridTemplateColumns.split(" ").length)).toBe(1);
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  }

  await page.setViewportSize({ width: 390, height: 844 });
  for (const route of workspaces) {
    await page.goto(route);
    await expect(page.locator(".design-layout").first()).toBeVisible();
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  }
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
