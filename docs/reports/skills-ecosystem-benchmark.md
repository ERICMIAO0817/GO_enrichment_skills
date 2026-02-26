# Skills 生态仓库对照（以 awesome-agent-skills 为入口）

日期：2026-02-26

## 1. 调研范围

从 `VoltAgent/awesome-agent-skills` 列表抽样了以下仓库进行结构与规范对照：

- `anthropics/skills`
- `openai/skills`
- `cloudflare/skills`
- `vercel-labs/agent-skills`
- `google-labs-code/stitch-skills`
- `huggingface/skills`
- `trailofbits/skills`
- `stripe/ai`
- `getsentry/skills`
- `microsoft/skills`

## 2. 观察到的共性模式

### A. Skill 本体模式（几乎所有仓库一致）

- 每个 skill 目录至少有 `SKILL.md`
- 常见附加目录：`scripts/`, `references/`, `assets/`, `agents/`
- `SKILL.md` 的最小 frontmatter 通常是：
  - `name`
  - `description`

### B. 发布与分发模式（高频）

- 根目录有清晰的 `README.md`（安装方式 + skill 列表）
- 有 `CONTRIBUTING.md`、`LICENSE`、`SECURITY.md`（至少其一）
- 多 agent 适配文件逐步成为标准：
  - `.claude-plugin/*`
  - `.cursor-plugin/*`
  - `AGENTS.md`（fallback 指令集）

### C. 可维护性模式（中高频）

- 自动生成 catalog（skill 表格 / AGENTS / plugin manifest）
- 通过脚本做发布工件再生成（例如 `publish.sh`）
- 强调“不要硬编码绝对路径”

## 3. 针对本仓库的落地原则

结合本仓库目标（GO 富集可靠性 skill）采用以下策略：

1. 保留单 skill 聚焦，不做过度 marketplace 复杂化。
2. 补齐发布面基础设施（README/CONTRIBUTING/SECURITY/LICENSE）。
3. 增加跨 agent 可发现文件（Claude/Cursor/AGENTS fallback）。
4. 用脚本生成/更新发布工件，减少手工漂移。
5. 保持 SKILL frontmatter 简洁，避免引入不必要字段。

## 4. 与调研仓库相比的取舍

- 未引入大型插件市场分层（如多 plugin marketplace），因为当前只有 1 个主技能。
- 未引入远程 MCP 清单（`.mcp.json`），因为当前 skill 不依赖自建远程 MCP endpoint。
- 先以“可发布 + 可维护 + 可移植”为第一阶段目标，后续再扩展多 skill 目录体系。

## 5. 下一步可扩展方向

- 若新增 3+ skills：引入自动生成 skill 目录索引与版本化发布脚本。
- 若要发布到 Claude/Cursor marketplace：补充更完整的 plugin metadata 与图标资源。
- 若要做质量门禁：增加 CI 对 `SKILL.md` frontmatter 和路径可移植性做 lint。
